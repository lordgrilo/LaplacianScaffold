import scipy
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import numpy as np
from scipy import linalg, matrix
import networkx as nx

def null(A, eps=1e-12):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)

def gs(X, row_vecs=False, norm = True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T



def IPR(v,p=4,verbose=False):
    n=float(np.linalg.norm(v))
    if n>0:
        if n!=1.0:
            if verbose==True:
                print 'Input vector with norm', n
            v=v/n
        return np.sum(map(lambda x: x**p,v));
    else:
        return IPR(np.ones(len(v)))


def kernel_base_IPR(x,nA):
    """
    the function expects as an input a vector with length equal to the 
    dimension of the kernel
    """
    return IPR(np.dot(nA,np.transpose(x)))

   
def orthogonality_bit(x,nA,new_base_matrix_list,mode='max'):
    s=[];
    for l in range(new_base_matrix_list.shape[1]):
        if np.linalg.norm(new_base_matrix_list[:,l])>0:
            s.append(np.dot(np.dot(nA,np.transpose(x))/np.linalg.norm(np.dot(nA,np.transpose(x))) , new_base_matrix_list[:,l]/np.linalg.norm(new_base_matrix_list[:,l])));
        else:
            s.append(0);
    if mode=='max':
        z = np.max(map(lambda x: x**2,s))
    if mode=='mean':
        z = np.mean(map(lambda x: x**2,s))
    return z

def find_optimal_base(nA,alpha,max_attempts=1, max_iter=100, pop_size=50, verbose=False):
    import copy
    strategies = ['best1bin', 'best1exp', 'rand1exp', 'randtobest1exp', 'best2exp', 'rand2exp', 'randtobest1bin', 'best2bin', 'rand2bin', 'rand1bin']

    attempt=0;
    if verbose==True:
        print 'Dimensions of null space', nA.shape
    base_matrix_list = nA   
    while attempt<max_attempts:
        attempt+=1;
        new_base_matrix_list=np.zeros_like(base_matrix_list);
        if verbose==True:
            print 'attempt', attempt
        #run over the kernel vectors 
        for k in range(nA.shape[1]):
            print 'Starting optimisation of vector ', k
            x0 = []
            for i in range(nA.shape[1]):
                x0.append([-1,1])
            x0=np.array(x0)
            inner_new_base_matrix_list=new_base_matrix_list[:,:k];
            def single_snapshot_trades(x):
                if inner_new_base_matrix_list.shape[1]>=1:
                    t = -alpha*kernel_base_IPR(x,nA) ## maximal IPR contributes to the minimazation 
                    r = (1-alpha)*orthogonality_bit(x,nA,inner_new_base_matrix_list);
                    return t+r
                else:
                    return -alpha*kernel_base_IPR(x,nA)
            res = differential_evolution(single_snapshot_trades, x0, strategy=strategies[attempt%len(strategies)] ,maxiter=max_iter,popsize=pop_size)#, T=0.1, niter=100,stepsize=.01)
            new_v = np.dot(nA,np.transpose(res['x']))/np.linalg.norm(np.dot(nA,np.transpose(res['x'])))
            new_base_matrix_list[:,k]=np.squeeze(new_v);

        if np.mean(map(lambda x: IPR(new_base_matrix_list[:,x]), range(new_base_matrix_list.shape[1])))>=np.mean(map(lambda x: IPR(base_matrix_list[:,x]), range(base_matrix_list.shape[1]))):
            print 'Updating base matrix with new one.'
            base_matrix_list = copy.copy(new_base_matrix_list);
            if verbose==True:
                print 'IPR report:', np.mean(map(lambda x: IPR(base_matrix_list[:,x]), range(new_base_matrix_list.shape[1]))), np.mean(map(lambda x: IPR(nA[:,x]), range(nA.shape[1])))
        else:
            print 'New base matrix is not satisfactory.'
            
    if verbose==True:
        print 'Final IPR report:', np.mean(map(lambda x: IPR(base_matrix_list[:,x]), range(new_base_matrix_list.shape[1]))), np.mean(map(lambda x: IPR(nA[:,x]), range(nA.shape[1])))


    return base_matrix_list



def recreate_scaffold(edge_weights,relabeled_simplex_ordered_base,dim):
    scaffold = nx.Graph()
    scaffold.add_nodes_from(map(lambda x: int(eval(x)[0]), relabeled_simplex_ordered_base[1].keys()))
    for edge_name in relabeled_simplex_ordered_base[dim]:
        edge =list(eval(edge_name))
        scaffold.add_edge(int(edge[0]), int(edge[1]),weight=edge_weights[relabeled_simplex_ordered_base[2][edge_name]])        
#        scaffold.add_edge(relabeled_simplex_ordered_base[1][str([edge[0]])], relabeled_simplex_ordered_base[1][str([edge[1]])],weight=edge_weights[relabeled_simplex_ordered_base[2][edge_name]])
    return scaffold;
