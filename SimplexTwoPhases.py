# Warning: Hàm chooseStarter bị sai nếu trong trường hợp cơ sở của phương án cực biên tìm được có chứa vector giả.
# Warning: Hàm giải thiếu trường hợp nghiệm bị suy biến.
# I don't have time to fix them so beware :).
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt
from tabulate import tabulate
@dataclass
class SimplexAlgorithmPackage:
    C:npt.NDArray[np.float16]
    # A:npt.NDArray[np.float16]
    B_Index:npt.NDArray[np.float16]
    # B:npt.NDArray[np.float16]
    # x:npt.NDArray[np.float16]
    Z:npt.NDArray[np.float16]
    Phi:npt.NDArray[np.float16]
    Delta:npt.NDArray[np.float16]
    f:float
    C_B:npt.NDArray[np.float16]
    X_B:npt.NDArray[np.float16]
    iter:int
def SimplexAlgorithm_check(A,b):
    (m,n)=np.shape(A)
    (h,k)=np.shape(b)
    if (n>=m and m==h and k==1):
        return
    raise TypeError(f"Matrix A and b missmatch. A must be of size mxn and b must be of size mx1, with n>=m.\n A is currently of size {np.shape(A)} and b is currently of size {np.shape(b)} ")
def SimplexAlgorithm_chooseStarter(C,A,b):
    (m,n)=np.shape(A)
    newC=np.zeros((1,n+3))
    newX=np.zeros((n+3,1))
    newA=np.ndarray.tolist(A)
    temp=[0.0]*m
    for i in range(0,3):
        newC[0,n+i]=1.0
        temp[i]=1.0
        newA[i]=newA[i]+temp
        temp[i]=0.0
        newX[n+i,0]=b[i,0]
    temp=SimplexAlgorithm_withStarter(newC,np.array(newA),b,newX,False)
    if(temp[1]!=0):
        raise Exception("Function doesn't have any solutions")
    return np.transpose(np.array([temp[0][0:n,0]]))
# return a tuple of (list of A_k,list of A_j) satisfying x_k=0,x_j!=0
def SimplexAlgorithm_Extreme_filterCorrespondingColumn(A:npt.NDArray[np.float16],x:npt.NDArray[np.float16])->tuple[npt.NDArray[np.float16],npt.NDArray[np.float16]]:
    neqB:list[npt.NDArray[np.float16]]=[]
    B:list[npt.NDArray[np.float16]]=[]
    for i in range(np.shape(A)[1]):
        if(x[i,0]>0):
            B.append(A[:,i])
        elif(x[i,0]==0):
            neqB.append(A[:,i])
        else:
            raise ValueError("x has a negative coordinate.")
    return (np.transpose(np.array(B)),np.transpose(np.array(neqB)))
# Only used in the first loop of the algorithm. Return Z_Bk
def SimplexAlgorithm_Column_findCoordinate(A_k:npt.NDArray[np.float16],B:npt.NDArray[np.float16])->npt.NDArray[np.float16]:
    return np.linalg.solve(B,A_k)
def SimplexAlgorithm_Column_findCoordinates(A,B):
    temp=[]
    for A_k in np.transpose(A):
        temp.append(SimplexAlgorithm_Column_findCoordinate(A_k,B))
    return np.transpose(np.array(temp))
def SimplexAlgorithm_Delta_cal(C,C_B,Z):
    """
    return an array of size 1xn with each entry has the value of delta_k
    """
    temp=[]
    for k in range(np.shape(C)[1]):
        temp.append(np.dot(C_B[0],Z[:,k])-C[0,k])
    return np.array([temp])
def SimplexAlgorithm_printState(Package:SimplexAlgorithmPackage):
    C_B=np.ndarray.tolist(Package.C_B)[0]
    B_Index=np.ndarray.tolist(Package.B_Index)[0]
    X_B=np.ndarray.tolist(np.transpose(Package.X_B))[0]
    SimplexBoard:list[list[str|np.float16|float]]=[["Col","C_B"]+C_B,["Base","B"]+[f"A{i+1}" for i in Package.B_Index[0]],["Sol","X"]+X_B] # type: ignore
    FinalLine:list[str|np.float16|float]=["",f"Board {Package.iter}",f"f={Package.f}"]+np.ndarray.tolist(Package.Delta)[0]+[""]
    for i in range(np.shape(Package.C)[1]):
        SimplexBoard.append([Package.C[0,i].astype(float),f"A{i+1}"]+np.ndarray.tolist(Package.Z[:,i]))
    SimplexBoard.append(["","Phi"]+np.ndarray.tolist(Package.Phi)[0])
    print(tabulate(np.ndarray.tolist(np.transpose(np.array(SimplexBoard)))+[[""]*len(FinalLine)]+[FinalLine]))
def SimplexAlgorithm_inner(Package:SimplexAlgorithmPackage,printOut:bool):
    n=np.shape(Package.Delta)[1]
    for k in range(n):
        if(Package.Delta[0,k]>0):
            phis:list[float]=[]
            minj=0
            for j in range(np.shape(Package.Z)[0]):
                if(Package.Z[j,k]>0):
                    phis.append(Package.X_B[j,0]/Package.Z[j,k])
                    if(phis[j]<phis[minj]):
                        minj=j
                else:
                    phis.append(float('inf'))
            if phis[minj]==float('inf'):
                raise Exception("Function doesn't have any solutions")
            
            Package.Phi=np.array([phis])
            if(printOut):
                SimplexAlgorithm_printState(Package)
            Package.C_B[0,minj]=Package.C[0,k]
            Package.B_Index[0,minj]=k
            Package.X_B[minj,0]=Package.X_B[minj,0]/Package.Z[minj,k]
            Package.Z[minj,:]=Package.Z[minj,:]/Package.Z[minj,k]
            for j in range(np.shape(Package.C_B)[1]):          
                if j==minj:
                    continue
                Package.X_B[j,0]=Package.X_B[j,0]-Package.Z[j,k]*Package.X_B[minj,0]
                Package.Z[j,:]=Package.Z[j,:]-Package.Z[j,k]*Package.Z[minj,:]
            Package.f=Package.f-Package.Delta[0,k]*Package.X_B[minj,0]
            Package.Delta[0,:]=Package.Delta[0,:]-Package.Delta[0,k]*Package.Z[minj,:]
            Package.iter=Package.iter+1
            # SimplexAlgorithm_printState(Package)
            # return
            return SimplexAlgorithm_inner(Package,printOut)
    temp=["" for i in range(np.shape(Package.Phi)[1])]
    Package.Phi=np.array([temp])
    if(printOut):
        SimplexAlgorithm_printState(Package)
    sol=[0.0 for i in range(n)]
    for i in range(np.shape(Package.X_B)[0]):
        sol[Package.B_Index[0,i]]=Package.X_B[i,0]
    return (np.transpose(np.array([sol])),Package.f)
         
                

def SimplexAlgorithm_withStarter(C,A,b,x,printOut=True):
    """
        C is a 1xn matrix\\
        A is a mxn matrix\\
        b is a mx1 matrix\\
        x is a nx1 matrix
    """
    SimplexAlgorithm_check(A,b)
    (B,neqB)=SimplexAlgorithm_Extreme_filterCorrespondingColumn(A,x)
    B_Index=SimplexAlgorithm_Extreme_filterCorrespondingColumn(np.array([[i for i in range(0,np.shape(A)[1])]]),x)[0]
    # print(C,x)

    Z=SimplexAlgorithm_Column_findCoordinates(A,B)
    (C_B,neqC)=SimplexAlgorithm_Extreme_filterCorrespondingColumn(C,x)
    Delta=SimplexAlgorithm_Delta_cal(C,C_B,Z)
    f=np.matmul(C,x)[0,0]
    X_B=np.transpose(SimplexAlgorithm_Extreme_filterCorrespondingColumn(np.transpose(x),x)[0])
    # Phi is first initialized as B_Index (only as a placeholder) 
    Package= SimplexAlgorithmPackage(C,B_Index,Z,B_Index,Delta,f,C_B,X_B,1)
    return SimplexAlgorithm_inner(Package,printOut)

def SimplexAlgorithm(C,A,b):
    return SimplexAlgorithm_withStarter(C,A,b,SimplexAlgorithm_chooseStarter(C,A,b))
C=np.array([[-2.0,2,1,3,0,1]])
A=np.array([[-1.0,2,2,0,1,0],[1,1,1,0,0,1],[2,0,2,1,0,0]])
b=np.transpose(np.array([[4.0,5,7]]))
# x=np.transpose(np.array([[0.0,0,0,7,4,5]]))
SimplexAlgorithm(C,A,b)
    
    