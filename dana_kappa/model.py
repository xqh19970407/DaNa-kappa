"""DANA additions to a separately installed ALiEGNN backbone.

Preserves the published checkpoint key layout, including dormant baseline
branches. The learned weights therefore load strictly without conversion.
"""
import math
import types

import dgl
import dgl.function as fn
import torch
from torch import nn


class NeighborhoodSlotAttention(nn.Module):
    def __init__(self, node_features=256, message_features=256):
        super().__init__()
        self.num_slots=4;self.slot_width=32;self.scale=self.slot_width**-.5
        self.message_key=nn.Linear(message_features,32,bias=False)
        self.message_value=nn.Linear(message_features,32,bias=False)
        self.node_query=nn.Linear(node_features,32,bias=False)
        self.message_query=nn.Linear(message_features,32,bias=False)
        self.slot_key=nn.Linear(32,32,bias=False)
        self.output=nn.Linear(32,message_features,bias=False)
        self.slot_vectors=nn.Parameter(torch.randn(4,32))
        self.alpha=nn.Parameter(torch.tensor(math.log(.1/.9)))

    def forward(self,graph,node_features,messages):
        if graph.num_edges()==0:return messages
        with graph.local_scope():
            values=self.message_value(messages)
            graph.edata['dana_key']=self.message_key(messages)
            query=self.node_query(node_features)
            logits=[]
            for s in range(4):
                graph.ndata['dana_query']=query+self.slot_vectors[s]
                graph.apply_edges(fn.e_dot_v('dana_key','dana_query','dana_score'))
                logits.append(graph.edata['dana_score']*self.scale)
            attention=dgl.ops.edge_softmax(graph,torch.cat(logits,dim=-1))
            slots=[]
            for s in range(4):
                graph.edata['dana_value']=values*attention[:,s:s+1]
                graph.update_all(fn.copy_e('dana_value','m'),fn.sum('m','dana_slot'))
                slots.append(graph.ndata['dana_slot'])
            states=torch.stack(slots,dim=1)
            keys=self.slot_key(states)
            graph.edata['dana_edge_query']=self.message_query(messages)
            logits=[]
            for s in range(4):
                graph.ndata['dana_slot_key']=keys[:,s,:]
                graph.apply_edges(fn.e_dot_v('dana_edge_query','dana_slot_key','dana_score'))
                logits.append(graph.edata['dana_score']*self.scale)
            attention=torch.softmax(torch.cat(logits,dim=-1),dim=-1)
            context=torch.zeros_like(values)
            for s in range(4):
                graph.ndata['dana_slot_value']=states[:,s,:]
                graph.edata['dana_read_weight']=attention[:,s:s+1]
                graph.apply_edges(fn.v_mul_e('dana_slot_value','dana_read_weight','dana_context'))
                context=context+graph.edata['dana_context']
            return messages+torch.sigmoid(self.alpha)*self.output(torch.nn.functional.silu(context-values))


def _aggregate(graph,messages):
    graph.edata['dana_m']=messages
    summaries=[]
    for reducer,field in [(fn.sum,'sum'),(fn.mean,'mean'),(fn.max,'max')]:
        graph.update_all(fn.copy_e('dana_m','m'),reducer('m','dana_'+field))
        summaries.append(graph.ndata['dana_'+field])
    graph.edata['dana_m2']=messages.square()
    graph.update_all(fn.copy_e('dana_m2','m'),fn.mean('m','dana_mean2'))
    variance=(graph.ndata['dana_mean2']-summaries[1].square()).clamp(min=0)
    summaries.append(torch.sqrt(variance+1e-6)-1e-3)
    present=(graph.in_degrees()>0).unsqueeze(-1)
    return torch.cat([torch.where(present,x,torch.zeros_like(x)) for x in summaries],dim=-1)


def _dana_forward(self,graph,node_feat,coord_feat,edge_feat=None):
    """Keep backbone geometry/edge updates, changing only node message reduction."""
    with graph.local_scope():
        graph.ndata['h']=node_feat;graph.ndata['x']=coord_feat
        if edge_feat is not None:graph.edata['a']=edge_feat
        graph.apply_edges(fn.u_sub_v('x','x','x_diff'))
        distance=torch.norm(graph.edata['x_diff'],dim=1,keepdim=True)
        graph.edata['x_vec']=graph.edata['x_diff']/(distance+1e-15)
        sh=self.sh_embedding(graph.edata['x_vec'],distance)
        graph.edata['sh_invariant']=self.readout(sh,self.sh_embedding)
        graph.apply_edges(self.message)
        raw=graph.edata['msg_h']
        corrected=self.neighborhood_attention(graph,node_feat,raw)
        node=self.node_network(torch.cat([node_feat,_aggregate(graph,corrected)],dim=-1))
        graph.update_all(fn.copy_e('msg_x','m'),fn.mean('m','coord_delta'))
        coord=coord_feat+graph.ndata['coord_delta']
        return (node,coord,edge_feat+raw) if edge_feat is not None else (node,coord)


def _layernorm(module):
    for name,child in list(module.named_children()):
        if isinstance(child,nn.BatchNorm1d):
            setattr(module,name,nn.LayerNorm(child.num_features,eps=1e-5,elementwise_affine=True))
        else:_layernorm(child)


def create_model(config):
    from aliegnn.models.alignn_egnn import ALIEGNN,ALIEGNNConfig
    model=ALIEGNN(ALIEGNNConfig(**config))
    _layernorm(model)
    convs=[block.node_update for block in model.GatedE3NNGCN_layers]+list(model.egcn_layers)
    for conv in convs:
        conv.node_network[0]=nn.Linear(conv.in_size+4*conv.hidden_size,conv.hidden_size)
        conv.neighborhood_attention=NeighborhoodSlotAttention(conv.in_size,conv.hidden_size)
        conv.forward=types.MethodType(_dana_forward,conv)
    return model


def loss(prediction,target):
    r=prediction.reshape(-1)-target.reshape(-1)
    return torch.mean(torch.hypot(r,r.new_tensor(.2))-.2)
