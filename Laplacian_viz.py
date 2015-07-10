import matplotlib.pyplot as plt 
import numpy as np
import networkx as nx


def draw_weighted_labeled_network(G,title='',pos=None):
    fig = plt.figure(figsize=(8,4))
    if pos==None:
        pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G,pos=pos);
    nx.draw_networkx_labels(G,pos=pos);
    nx.draw_networkx_edges(G,pos=pos);
    nx.draw_networkx_edge_labels(G,pos,edge_labels=nx.get_edge_attributes(G,'weight'))
    plt.title(title);
    plt.tight_layout()
    plt.show()


def draw_kernel_eigenvector_reduced(scaffold,pos,relabeled_simplex_ordered_base,dim,kernel_eigenvector):
    g = nx.Graph();
#    g.add_nodes_from(scaffold.nodes());
#    g.add_edges_from(scaffold.edges());
    w_dict = {}
    nodes = []
    for edge_name in relabeled_simplex_ordered_base[dim]:
        if relabeled_simplex_ordered_base[dim][edge_name] in kernel_eigenvector.index:
            edge = list(eval(edge_name))
            w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
            w_dict[(int(edge[0]),int(edge[1]))] = w;
            nodes.append(int(edge[0]))
            nodes.append(int(edge[1]))
    g.add_nodes_from(nodes)
    g.add_edges_from(w_dict.keys());
    nx.draw_networkx(g,pos,node_color='r',node_size=300,alpha=0.8)
    nx.draw_networkx_edges(g,pos,width=1.0,alpha=0.5);    

    for edge_name in relabeled_simplex_ordered_base[dim]:
        if relabeled_simplex_ordered_base[dim][edge_name] in kernel_eigenvector.index:
            edge = list(eval(edge_name))
            w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
            nx.draw_networkx_edges(g,pos,edgelist=[(int(edge[0]),int(edge[1]))],width=50*w,alpha=0.5,edge_color='r')
    return



def draw_kernel_eigenvector(scaffold,pos,relabeled_simplex_ordered_base,dim,kernel_eigenvector,factor_w=100,nodesize=100):
    g = nx.Graph();
    g.add_nodes_from(scaffold.nodes());
    g.add_edges_from(scaffold.edges());
    nx.draw_networkx(g,pos,node_color='r',node_size=nodesize,alpha=0.8)
    #nx.draw_networkx_edges(g,pos,width=1.0,alpha=0.5);    

    for edge_name in relabeled_simplex_ordered_base[dim]:
        if relabeled_simplex_ordered_base[dim][edge_name] in kernel_eigenvector.index:
            edge =list(eval(edge_name))
            w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
            nx.draw_networkx_edges(g,pos,edgelist=[(int(edge[0]),int(edge[1]))],width=factor_w*w,alpha=0.5,edge_color='r')
    return 


def draw_aggregated_kernel_eigenvector(scaffold,pos,relabeled_simplex_ordered_base,dim,kernel_eigenvector,factor_w=100):
    g = nx.Graph();
    g.add_nodes_from(scaffold.nodes());
    g.add_edges_from(scaffold.edges());
    nx.draw_networkx(g,pos,node_color='r',node_size=300,alpha=0.8)
    nx.draw_networkx_edges(g,pos,width=1.0,alpha=0.5);    

    for edge_name in relabeled_simplex_ordered_base[dim]:
        edge =list(eval(edge_name))
        w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
        nx.draw_networkx_edges(g,pos,edgelist=[(int(edge[0]),int(edge[1]))],width=factor_w*w,alpha=0.5,edge_color='b')
    return
