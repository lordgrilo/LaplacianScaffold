
# coding: utf-8

# In[1]:


import pickle
import os, sys
import pandas as pd
from numpy import *
import numpy as np

# FILTRATION GIVEN BY THE USER, MUST BE A DICTIONARY WITH KEYS ALL THE SIMPLICES IN THE SIMPLICIAL COMPLEX AND AS VALUE THE STEP OF THE FILTRATION WHEN THEY ARE CREATED (AND A WEIGHT NOT NECESSARY)
#test_filtration=sys.argv[1];
#output_bm_file=sys.argv[2];
#fil=pickle.load(test_filtration);



# FUNCTION CODIMENSION: GIVEN 2 SIMPLICES RETURNS 1 AND -1 IF THAT IS THAT IS THEIR CODIMENSION AND 0 OTHERWISE
def codimension(simplex1,simplex2,verbose=False):

    if verbose==True:
         ( a, len(a), b, len(b))
    if (simplex1 < simplex2) and (len(simplex2-simplex1)==1):
        return 1;
    if (simplex2 < simplex1) and (len(simplex1-simplex2)==1):
        return -1;
    return 0;

# FUNCTION: INVERTS THE DICTIONARY GIVEN SO THAT THE KEYS ARE COUPLES WHERE THE FIRST NUMBER IS THE FILTRATION STEP AND THE SECOND IS THE DIMENSION+1 OF THE SIMPLICES IN THE CORRESPONDING VALUE

def invert_filtration_dictionary(fil):
    import ast
    inv_dict={};
    
    for c in fil: # TAKES A SIMPLEX IN THE KEYS OF FILTRATION
        try:
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c))); # int(fil[c][0]) IS THE FILTRATION STEP WHERE c IS ADDED, AND len(set(ast.literal_eval(c)))) IS HIS LENGHT. IT LOCATES THIS KEY IN THE NEW FILTRATION DICTIONARY AND APPENDS THE SIMPLEX c TO THE VALUE.
        except: # IF THE KEY DOESN'T ALREADY EXIST
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))]=[]; # CREATES IT
            inv_dict[(int(fil[c][0]),len(set(ast.literal_eval(c))))].append(set(ast.literal_eval(c))); # AND APPENDS THE SIMPLEX c TO THE VALUE.
    for l in inv_dict:
        inv_dict[l]=sorted(inv_dict[l], key=lambda x: len(x)); # SORTS THE KEYS FOR LENGHT (USEFUL?)
        
    return inv_dict;

# FUNCTION SPARSE_BOUNDARY_MATRIX_ZERO: RETURNS THE kth-BOUNDARY MATRIX FOR THE SIMPLICIAL COMPLEX AT STEP c IN THE FILTRATION WHERE THE FIRST k-1-SIMPLEX IS ADDED TO THE FILTRATION (YOU DON'T NEED TO ADD THE PREVIOUS kth-BOUNDARY MATRIX)
def sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False): 
    
    ordered_simplex_list=[]; # THE LIST WITH ALL THE SIMPLICES IT NEEDS TO CHECK TO CREATE THE MATRIX
    Ord=[] # LIST OF THE (k-1)-SIMPLICES THAT ARE IN THE SIMPLICIAL COMPLEX TO PASS TO THE NEXT STEP IN THE FILTRATION
    try:
        ordered_simplex_list.extend(inv_dict[(c,k)]); # ADDS THE (k-1)-SIMPLICES ADDED AT STEP c TO THE LIST
        R=len(ordered_simplex_list); # NUMBER OF ROWS IN THE MATRIX
        Ord=list(ordered_simplex_list) # CREATES THE LIST TO BE RETURNED
	#print 'here i am born %d'% c,Ord
    except KeyError:
	#print 'ord sono al primo try',Ord
        return matrix([]),Ord # IF THERE ARE NO (k-1)-SIMPLICES TO ADD, SINCE THERE WHERE NONE IN THE PREVIOUS STEPS EITHER IT RETURN AN EMPTY MATRIX AND AN EMPTY LIST
    try:
        ordered_simplex_list.extend(inv_dict[(c,k+1)]);  # ADDS THE k-SIMPLICES ADDED AT STEP c TO THE LIST
    except KeyError:
	#print 'ord sono al secondo try',Ord
        return matrix([]),Ord # IF THERE ARE NO k-SIMPLICES TO ADD, THE MATRIX HAS NO COLUMNS SO IT RETURNS AN EMPTY MATRIX AND THE LIST OF (k-1)-SIMPLICES ADDED AT THIS STEP
    C=len(ordered_simplex_list)-R; # THE NUMBER OF COLUMNS IN THE MATRIX
    ordered_simplex_series=pd.Series(ordered_simplex_list,index=range(len(ordered_simplex_list))); # CREATES A PANDAS.SERIES TO USE AS REFERENCE GUIDE WHEN CREATING THE MATRIX
    del ordered_simplex_list;
    bm=zeros((R,C)); # CREATES A ZERO MATRIX OF THE RIGHT DIMENSIONS
    for i in ordered_simplex_series.index: # TAKES A NUMER IN THE RANGE
        if len(ordered_simplex_series[i])==k: # IF THE CORRESPONDING SIMPLEX HAS DIMENSION k
            for j in ordered_simplex_series.index[i:]: # FOR THE SIMPLICES AFTER THAT CHECKS THE CODIMENSION
                cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
                if cod==1:
		    piu=list(ordered_simplex_series[j]-ordered_simplex_series[i])
		    s=sorted(ordered_simplex_series[j])
		    esp=s.index(piu[0])
		    if esp%2==0:
		    	bm[i,j-R]=1;
		    else:
			bm[i,j-R]=-1;
                    if verbose==True:
                        print (i,j, cod,ordered_simplex_series[i], ordered_simplex_series[j]);
    bm= matrix(bm);
    #print 'sum col',[sum([x]) for x in bm.columns()]
    #print 'ord sono alla fine',Ord
    return matrix(bm),Ord

def sparse_boundary_matrix(inv_dict,c,k,deltak=matrix([]),ordered_ksimplex_list=[]):#where deltak is the kth boundary matrix of c-1
    if c==0:
        return sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False)
    elif not ordered_ksimplex_list:
	print 'I have found that there are no %d simplices up to this step'%(k-1),ordered_ksimplex_list;
	return sparse_boundary_matrix_zero(inv_dict,c,k,verbose=False)
    ordered_simplex_list=[];
    new=False;
    Ord=[]
    Ord=list(ordered_ksimplex_list);
    try:
        ordered_simplex_list.extend(inv_dict[(c,k)]);
        R=len(ordered_simplex_list);
        Ord.extend(ordered_simplex_list)
        if deltak.shape!=(1,0):
            D=matrix(vstack([deltak, matrix(zeros((R,shape(deltak)[1])))]))
	    row_sum_D=[sum([x]) for x in D]
	else:
	    row_sum_D=[0]*(len(Ord))
    except KeyError:
        R=0;
        print 'no new %d simplex'%(k-1) 
        D=deltak;
	row_sum_D=[sum([x]) for x in D]

    try:
        ordered_simplex_list.extend(inv_dict[(c,k+1)]);
        if deltak.shape==(1,0):
            new=True;
    except KeyError:
        if deltak.shape==(1,0):
            return matrix([]),Ord
        else:
            return matrix(D),Ord
    C=len(ordered_simplex_list)-R;
    ordered_ksimplex_list.extend(ordered_simplex_list)
    ordered_simplex_series=pd.Series(ordered_ksimplex_list,index=range(len(ordered_ksimplex_list)));
    del ordered_ksimplex_list;
    del ordered_simplex_list;
    bm=zeros((len(Ord),C));
    for i in ordered_simplex_series.index:
        if len(ordered_simplex_series[i])==k:
            for j in ordered_simplex_series.index[i:]:
                cod=codimension(ordered_simplex_series[i], ordered_simplex_series[j]);
                if cod==1:
		    piu=list(ordered_simplex_series[j]-ordered_simplex_series[i])
		    s=sorted(ordered_simplex_series[j])
		    esp=s.index(piu[0])
		    if esp%2==0:
		    	bm[i,j-len(Ord)]=1;
		    else:
			bm[i,j-len(Ord)]=-1;
		
    if new:
	bm= matrix(bm);
        BM=bm;
        del bm
    else:
        BM=hstack([D,bm]);
	BM= matrix(BM);
        del bm,D
    #print [sum([x]) for x in BM.columns()]
    return matrix(BM),Ord;
        
    
def Laplacian(inv_fil,c,k,deltak=matrix([]),Ord_k=[],deltak1=matrix([]),Ord_k1=[],save_boundary=True,verbose=False):
    if k==0:
        Dk1,Ordk1=sparse_boundary_matrix(inv_fil,c,k+1,deltak1,Ord_k1)
        print (c,k+1)
	if Dk1:
	   Dk1=matrix(Dk1); 
	   L=Dk1*(Dk1.transpose())
	   return L,matrix([]),Dk1,[],Ordk1
	else:
	   return matrix([]),matrix([]),Dk1,[],Dk1
    Dk,Ordk=sparse_boundary_matrix(inv_fil,c,k,deltak,Ord_k)
    print (c,k)
    Dk1,Ordk1=sparse_boundary_matrix(inv_fil,c,k+1,deltak1,Ord_k1)
    print (c,k+1) 
    if Dk1.shape==(1,0):
        if Dk.shape==(1,0):
            return matrix([]),Dk,Dk1,Ordk,Ordk1
        else:
	    if verbose:
	    	print 'dk\n', Dk;
	    Dk=matrix(Dk)
            L=(Dk.transpose())*Dk
    else:
	Dk=matrix(Dk)
	Dk1=matrix(Dk1)
	if verbose:
	    print 'dk\n', Dk;
	    print 'dk1\n', Dk1;
        L=(Dk.transpose())*Dk+Dk1*(Dk1.transpose())
    if save_boundary==True:
        return L,Dk,Dk1,Ordk,Ordk1
    else:
        return L

def Proiettore(L):
    LT=L.transpose()
    PP=LT*L
    invPP=(PP).inverse()
    P=((L*invPP)*LT)
    return P




def right_kernel_space(L):
    if L==0:
        return Matrix(ZZ,[])
    u, s, vh = scipy.linalg.svd(L)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    if null_space.any():
        return  Matrix(ZZ,null_space.transpose())
    else:
        return Matrix(ZZ,[])

def column_space(L):
    #print('COLUMN SPACE - I am using toll:', eps)
    if L==0:
        return Matrix(ZZ,[])
    u, s, vh = scipy.linalg.svd(L)
    column_mask = (s >= eps)
    column = scipy.compress(column_mask, u, axis=1)
    return  Matrix(ZZ,column)

def H(LambdaD, BD, BBD,verbose=False):
    k=[1,1,1]
    if verbose==True:
	print parent(LambdaD).dims(), parent(BD).dims(), parent(BBD).dims();
    if not LambdaD:
        k[0]=0
    if not BD:
        k[1]=0
    if not BBD:
        k[2]=0
    if verbose==True:
	print k;
    if k==[1,1,1]:
        M=hstack([hstack([LambdaD,BD]),BBD])
    elif k==[1,0,1]:
        M=hstack([LambdaD,BBD])
    elif k==[1,1,0]:
        M=hstack([LambdaD,BD])
    elif k==[1,0,0]:
        return LambdaD
    elif k==[0,1,1]:
        M=hstack([BD,BBD])
    elif k==[0,1,0]:
        return BD
    elif k==[0,0,1]:
        return BBD
    return M



def Persistent_Homology_maps(k,verbose=False):
    from numpy import zeros
    n=sorted(inv_fil.keys())[-1][0]
    homCD={}
    LapC,D1C,D2C,O1C,O2C=Laplacian(inv_fil,0,k)
    LambdaC=( Matrix(ZZ,LapC )).kernel()
    LambdaC=LambdaC.basis_matrix()
    LambdaC=(Matrix(ZZ,LambdaC)).transpose()
    del LapC
    for c in range(1,n+1):
        print c 
        LapD,D1D,D2D,O1D,O2D=Laplacian(inv_fil,c,k,D1C,O1C,D2C,O2C)
        if [LapD,D1D,D2D]!=[0,0,0]:
            if [D1C,D2C]==[0,0]:
                LambdaD=(Matrix(ZZ,LapD)).kernel()
		LambdaD=LambdaD.basis_matrix()
                LambdaC=(Matrix(ZZ,LambdaD)).transpose()
		print 'LambdaD con tutti 0',LambdaD;
                D1C=D1D
                if D2D!=0:
                    D2C=D2D
            else:
                LambdaD=(Matrix(ZZ,LapD)).kernel()
		LambdaD=(Matrix(ZZ,LambdaD.basis_matrix())).transpose()
		if verbose:
		    print 'LambdaD',LambdaD;
		PD=Proiettore(LambdaD)

		BD=(Matrix(ZZ,D2D)).column_space()
		BD=(BD.basis_matrix()).transpose()
		if verbose:
		    print 'Ord1',O1D,'\n Ord2',O2D;
		    print D2D
		    print D1D
                fuffa=(D1D).transpose()
                BBD=(Matrix(ZZ,(fuffa))).column_space()
		BBD=(BBD.basis_matrix()).transpose()

		if D1C:
		    lentC=shape(D1C)[1];
		else:
		    lentC=shape(D2C)[0]
		if D1D:
		    lentD=shape(D1D)[1];
		else:
		    lentD=shape(D2D)[0]
                IDC=eye(lentC);
                if lentD>lentC:
		   r=lentD-lentC
                   ZERO=zeros((r,lentC));
                   F1C=Matrix(ZZ,vstack([IDC,ZERO]));
                else:
                   F1C=Matrix(ZZ,IDC);
		    
                HD=H(LambdaD,BD,BBD)
		HD=Matrix(ZZ,HD)
		if  HD.is_square():
		    print 'HD is square';
		    HD=Matrix(ZZ,HD).inverse()
		else:
                    print 'HD=LambdaD|BD|BBD',(LambdaD.nrows(),LambdaD.ncols()),(BD.nrows(),BD.ncols()),(BBD.nrows(),BBD.ncols());
		    print BD
                    raise ValueError ('ERROR HD NOT SQUARE')

		if not LambdaC:
		    homCD[c-1]=Matrix(ZZ,[])
		    LambdaC=LambdaD
                    D1C=D1D
               	    if D2D!=0:
                    	D2C=D2D
		else:
		    if verbose:
		    	print 'HD \n',HD,'\n PD \n',PD,'\n F1C\n', F1C,'\n LambdaC \n',LambdaC;
			
                    HOM=HD*(PD*(F1C*LambdaC))
		    homCD[c-1]=HOM[:LambdaD.ncols()][:]
                    LambdaC=LambdaD
                    D1C=D1D
                    if D2D!=0:
                  	  D2C=D2D
        O1C=O1D#4
        O2C=O2D#
    if D1D:
	lentD=shape(D1D)[1];
    else:
	lentD=shape(D2D)[0]
    IDD=Matrix.identity(lentD);
    try:
    	HOM=HD*(PD*(IDD*LambdaD))
    except UnboundLocalError: #local variable 'HD' referenced before assignment
        del D1C,D2C,O1C,O2C,LapD,D1D,D2D,O1D,O2D,LambdaC
	return homCD
    if HOM:
	homCD[c]=HOM[:LambdaD.ncols()][:]
    #print 'THIS IS THE LAST STEP\n HD\n',HD,'\n PD\n',PD,'\n IDD\n',IDD,'\n LambdaD \n',LambdaD
    del D1C,D2C,O1C,O2C,LapD,D1D,D2D,O1D,O2D,LambdaC
    print 'Done.' 
    return homCD


import scipy 
from scipy.sparse.linalg import svds

def null(A, eps=1e-12,sparse=True):
    if sparse == True:
        X = scipy.sparse.csc_matrix(A)
        n=X.shape[1]
        u, s, vh = svds(X, n-1, which='SM')
    else:
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


def calculate_nullspaces(Laplacian_dict):
    nullspaces={}
    import time
    import numpy as np
    for l,e in enumerate(sorted(Laplacian_dict.keys())):
        if l>=1:
            internal_now=time.time()
            L=Laplacian_dict[e];
            L_old=Laplacian_dict[Laplacian_dict.keys()[l-1]];

            if (L.shape==L_old.shape):#  and (L==L_old).all(): unnecessary condition
                print 'No changes. Skipping step:' , e
                #nullspaces[e] = nullspaces[Laplacian_dict.keys()[l-1]]
            else:
                nA=null(Laplacian_dict[e])
                nullspaces[e] = gs(nA);
                print 'Step:', e, ' Dimension nullspace:', Laplacian_dict[e].shape[1], ' elapsed time ', time.time()-internal_now
        else:
            internal_now=time.time()
            nA=null(Laplacian_dict[e])
            nullspaces[e]=gs(nA);
            print 'Step:', e, ' Dimension nullspace:', Laplacian_dict[e].shape[1], ' elapsed time ', time.time()-internal_now
    return nullspaces;



def calculate_optimal_basis(inv_fil,k,null_spaces,relabeled_simplex_ordered_base,alpha=0.5,max_attempts=10,verbose=False):
    import copy, time
    import Laplacian_basis_optimization as lbo;

    dimensions=list(set(map(lambda x: x[1], inv_fil.keys())))
    steps=sorted(list(set(map(lambda x: x[0], inv_fil.keys()))))
    simplex_ordered_at_step={}
    simplex_ordered_at_step[0] = list(map(list,inv_fil[(0,k+1)]));
    optimal_bases={}

    # initialization of the first basis, for use later
    s=steps[0]
    print 0,s
    print 'Dimension of null space', null_spaces[s].shape

    if null_spaces[s].shape[1]>1:
        optimal_bases[s] = (lbo.find_optimal_base(null_spaces[s],alpha,max_attempts,verbose), map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
    elif null_spaces[s].shape[1]==1:
        optimal_bases[s] = (null_spaces[s]/np.linalg.norm(null_spaces[s]), map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
    else:
        optimal_bases[s] = ([], map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
    print '\n'
    
    ## calculation of all other steps, starting from step 1 
    for i,s in enumerate(steps[1:]):
        print i,s
        simplex_ordered_at_step[s] = copy.copy(simplex_ordered_at_step[steps[i]]);
        simplex_ordered_at_step[s].extend(list(map(list,inv_fil[s,k+1])))
        now=time.time()
        if null_spaces[s].shape[1]>1:
            optimal_bases[s] = (lbo.find_optimal_base(null_spaces[s],alpha,max_attempts,verbose), map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
        else:
            if null_spaces[s].shape[1]==1:
                optimal_bases[s] = (null_spaces[s]/np.linalg.norm(null_spaces[s]), map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
            else:
                optimal_bases[s] = ([], map(lambda x: relabeled_simplex_ordered_base[k+1][str(x)], list(simplex_ordered_at_step[s])));
        print '\n'

    return optimal_bases


def eigen_weights_creator(optimal_bases,mode):
    indexed_bases_dict={}
    for e in sorted(optimal_bases.keys()):
        if len(optimal_bases[e][0])>0:
            indexed_bases=pd.DataFrame(np.array(optimal_bases[e][0]))
            print e, pd.DataFrame(optimal_bases[e][0]).shape, len(optimal_bases[e][1])
            if indexed_bases.shape[1]>0:
                if mode=='quantum':
                    # quantum-wave-like version
                    indexed_bases_dict[e]=pd.DataFrame(np.power(optimal_bases[e][0].T,2), columns=optimal_bases[e][1])
                elif mode=='L1':
                    # L1-norm-like version
                    indexed_bases_dict[e]=pd.DataFrame(np.abs(optimal_bases[e][0]), columns=optimal_bases[e][1])
                else:
                    print 'Mode unspecified or invalid.'
    return indexed_bases_dict;



def recreate_scaffold(edge_weights,relabeled_simplex_ordered_base,dim):
    import networkx as nx;
    scaffold = nx.Graph()
    scaffold.add_nodes_from(map(lambda x: int(eval(x)[0]), relabeled_simplex_ordered_base[1].keys()))
    for edge_name in relabeled_simplex_ordered_base[dim+1]:
        edge =list(eval(edge_name))
        scaffold.add_edge(int(edge[0]), int(edge[1]),weight=edge_weights[relabeled_simplex_ordered_base[dim+1][edge_name]])        
    edges = scaffold.edges(data=True);
    for edge in edges:
        if edge[2]['weight']<=0:
            scaffold.remove_edge(edge[0],edge[1])

    return scaffold;


def Laplacian_scaffold(g,k=1,alpha=0.5, max_attempts=5, mode='quantum',verbose=False,filtration='standard'):
    import Holes as ho
    import time
    now=time.time()
    print 'Computing standard weight rank filtration.'
    if filtration=='standard':
        fil = ho.filtrations.standard_weight_clique_rank_filtration(g);
    elif filtration=='dense':
        fil = ho.filtrations.dense_graph_weight_clique_rank_filtration(g,k);
    elif filtration=='limited':
        fil = ho.filtration.limited_weight_clique_rank_filtration(g,k+2);

    print 'Completed.', time.time()-now;
    print 'Inverting filtration for Laplacian calculation.'
    now=time.time()
    inv_fil=invert_filtration_dictionary(fil)
    print 'Complete.', time.time()-now;

    print 'Sorting simplices in order of appearance and dimension.'
    now=time.time()
    simplex_ordered_base={}
    dimensions=list(set(map(lambda x: x[1], inv_fil.keys())))
    for d in sorted(dimensions):
        simplex_ordered_base[d]=[]
        for key in sorted(inv_fil.keys()):
            if key[1]==d:
                simplex_ordered_base[d].extend(map(list,inv_fil[key]))

    relabeled_simplex_ordered_base={}
    for d in dimensions:
        relabeled_simplex_ordered_base[d]={}
        for i,el in enumerate(simplex_ordered_base[d]):
            relabeled_simplex_ordered_base[d][str(el)]=i;
    print 'Done. Elapsed time: ', time.time()-now;  

    print 'Starting construction of Laplacians.'
    now=time.time();
    Laplacian_dict={};
    Laplacian_dict[0],D1C,D2C,O1C,O2C=Laplacian(inv_fil,0,k);
    for c in range(1,len(inv_fil)):
        print c 
        Laplacian_dict[c],D1C,D2C,O1C,O2C=Laplacian(inv_fil,c,k,D1C,O1C,D2C,O2C);
    print 'Done. Elapsed time: ', time.time()-now;  

    print 'Obtaining Laplacian nullspace basis.'
    now=time.time();
    nullspaces = calculate_nullspaces(Laplacian_dict);
    print 'Done. Elapsed time: ', time.time()-now;  

    print 'Starting optimization of nullspaces basis for Laplacian scaffold.'
    now = time.time()
    optimal_bases = calculate_optimal_basis(inv_fil,k,nullspaces,relabeled_simplex_ordered_base,alpha,max_attempts,verbose);
    print 'Done. Elapsed time: ', time.time()-now;  

    print 'Creating positive valued vectors for scaffold creations';
    print 'Chosen mode: ', mode;
    now = time.time()
    indexed_bases_dict = eigen_weights_creator(optimal_bases,mode);
    print 'Done. Elapsed time: ', time.time()-now;  

    print 'Weighting single snapshot contributions by slice width.'
    now = time.time()
    laplacian_weights = {}
    original_weights=sorted(list(set(nx.get_edge_attributes(g,'weight').values())),reverse=True)
    original_weights.append(0)
    w=np.diff(original_weights[::-1])[::-1];
    for i,e in enumerate(optimal_bases):
        laplacian_weights[e]=w[i];
    deformed_stepwise_vectors = pd.DataFrame(columns=indexed_bases_dict.keys(), index=(relabeled_simplex_ordered_base[2].values()));
    for e in indexed_bases_dict:    
        deformed_stepwise_vectors[e] = laplacian_weights[e]*indexed_bases_dict[e].sum()
    deformed_stepwise_vectors.fillna(0);

    ## final scaffold edge weights that can be directly used. 
    edge_weights = deformed_stepwise_vectors.sum(axis=1);
    print 'Done. Elapsed time: ', time.time()-now;  
    
    print 'Creating Laplacian scaffold.'
    now = time.time()
    scaffold = recreate_scaffold(edge_weights, relabeled_simplex_ordered_base,k)
    print 'Done. Elapsed time: ', time.time()-now;  
    return scaffold;


import matplotlib.pyplot as plt;
import networkx as nx

def draw_kernel_eigenvector(scaffold,pos,relabeled_simplex_ordered_base,dim,kernel_eigenvector):
    import networkx as nx
    g = nx.Graph();
    g.add_nodes_from(scaffold.nodes());
    g.add_edges_from(scaffold.edges());
    nx.draw_networkx(g,pos,node_color='r',node_size=300,alpha=0.8)
    nx.draw_networkx_edges(g,pos,width=1.0,alpha=0.5);    

    for edge_name in relabeled_simplex_ordered_base[dim]:
        if relabeled_simplex_ordered_base[dim][edge_name] in kernel_eigenvector.index:
            edge =list(eval(edge_name))
            w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
            nx.draw_networkx_edges(g,pos,edgelist=[(int(edge[0]),int(edge[1]))],width=50*w,alpha=0.5,edge_color='r')
    plt.show()
    return

def draw_aggregated_kernel_eigenvector(scaffold,pos,relabeled_simplex_ordered_base,dim,kernel_eigenvector):
    g = nx.Graph();
    g.add_nodes_from(scaffold.nodes());
    g.add_edges_from(scaffold.edges());
    nx.draw_networkx(g,pos,node_color='r',node_size=300,alpha=0.8)
    nx.draw_networkx_edges(g,pos,width=1.0,alpha=0.5);    
    for edge_name in relabeled_simplex_ordered_base[dim]:
        edge =list(eval(edge_name))
        w = kernel_eigenvector[relabeled_simplex_ordered_base[dim][edge_name]];
        nx.draw_networkx_edges(g,pos,edgelist=[(int(edge[0]),int(edge[1]))],width=10*w,alpha=0.5,edge_color='b')

    plt.show()
    return








