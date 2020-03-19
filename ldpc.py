import numpy as np
# Define binary calculation in rows
def row_cal(a,b):
    if len(a) != len(b):
        print('Input a and b with the same length')
        return
    list = []
    for i in range(len(a)):
        list.append(int(a[i])^int(b[i]))
    return np.array(list)    

### Encoding function ###
def Encode(H):
### Producing reduced row echolen matrix H_hat ###
    A = H.copy()
    M, N = A.shape
    i = 0
    row_exchanges = np.arange(M)
    for j in range(N):
        pivot = np.argmax(A[i:M, j]) + i
        val = A[pivot, j]
        if val == 0: #skip if the column is all zeros
            A[i:M, j] = np.zeros(M-i)
        else:
            if pivot != i: #whether to exchange rows?
                A[[pivot, i], j:N] = A[[i, pivot], j:N]
                #record changes of rows
                row_exchanges[[pivot,i]] = row_exchanges[[i,pivot]] 
            v = A[i, j:N] #ready to use for row calculation
            if i > 0: #if current row is not the first row
                A_above = A[np.arange(i),j:N] ## rows above the current row
                 #loop according to the first element of each row
                for idx,val in enumerate(A_above[:,0]):
                    if val == 0: 
                        A_above[idx] = A_above[idx] #if 0, then do calculation
                    else: 
                        A_above[idx] = row_cal(v, A_above[idx]) 
                A[np.arange(i),j:N]= A_above #replace rows above with calculated ones
            if i < M-1: #if current row is not the last row
                A_below = A[np.arange(i+1,M),j:N] #rows below the current row
                for idx,val in enumerate(A_below[:,0]): #loop according to the first element of each row
                    if val == 0: 
                        A_below[idx] = A_below[idx] #if 0, then do calculation
                    else: 
                        #if 1, then do calculation
                        A_below[idx] = row_cal(v, A_below[idx]) 
                #replace rows below with calculated ones        
                A[np.arange(i+1,M),j:N]= A_below 
            i += 1 #row+1
        if i == M: #finish all rows
            break
## Find I matrix and P matrix for G matrix ##
    row_idx = np.where(np.sum(A,axis=1)>0)[0]
    col_idx = np.where(np.sum(A,axis=0)>1)[0]
    P = A[row_idx,:]
    P = P[:,col_idx]
    B = np.vstack((np.eye(len(col_idx)), P))
    return A, col_idx, B
    
H_t = np.array([[1,1,1,1,0,0],[0,0,1,1,0,1],[1,0,0,1,1,0]])
H_hat_t, col_idx_t, G_t = Encode(H_t)
print("matrix H_hat_t is :\n",H_hat_t)
print("colums of H_hat to generate G:\n",col_idx_t+1)
print("matrix G_t is:\n",G_t)

### Define Deode function ###
def Decode(H_hat, y, p, max_iter=20, display=False):
    M, N = np.shape(H_hat)
    ms, ns = np.where(H_hat==1)
    syndrome = (H_hat @ y) % 2
    ### Initialization ###
    p1 = np.zeros(N)
    message_r = np.zeros(np.shape(H_hat))
    message_q = np.zeros(np.shape(H_hat))
    for i in range(N):
        if y[i]==1:
            message_q[ms[ns==i],i] = 1-p
            p1[i] = 1-p
        else:
            message_q[ms[ns==i],i] = p
            p1[i] = p
    post_hat = np.ones(np.shape(y))
    x_hat = np.ones(np.shape(y))>0.5
    status = -1
    if display:print("initial message r is:\n",message_r,
                     "\ninitial message q is:\n",message_q,
                     "\ninitial prior p1 is :\n",p1)
                     
    ### Message Update ###
    for i in range(max_iter):
    ### Horizontal: Message from factors to variables ###
        for m in range(M):
            n_idx = ns[ms==m]
            for n in n_idx:
                var_idx = n_idx[n_idx!=n]
                delta_q = (1-message_q[m,var_idx]) - message_q[m,var_idx]
                delta_r = np.prod(delta_q)
                message_r[m,n] = 0.5*(1-delta_r)
        if display: print("\nmessage_r of update %d is :\n"%(i+1), message_r)
        ### Vertical: Message from variables to factors ###
        for n in range(N):
            m_idx = ms[ns==n]
            for m in m_idx:
                fac_idx = m_idx[m_idx!=m]
                q1 = p1[n]*np.prod(message_r[fac_idx,n])
                q0 = (1-p1[n])*np.prod(np.ones(np.shape
                                               (message_r[fac_idx,n]))
                                       -message_r[fac_idx,n]) 
                message_q[m,n] = q1/(q1+q0)
                post1 = p1[n]*np.prod(message_r[m_idx,n])
                post0 = (1-p1[n])*np.prod(np.ones(np.shape(message_r[m_idx,n]))
                                          -message_r[m_idx,n])
                post_hat[n] = post1/(post1+post0)
                if post_hat[n]>0.5: x_hat[n] = 1
                else: x_hat[n] = 0
            if display:
                print("message_q of update %d is :\n"%(i+1),message_q)
                print("post_hat of update %d is :\n"%(i+1),post_hat)
                print("x_hat of update %d is :\n"%(i+1),x_hat)
                ### Stopping ###
            if np.all((H_hat@x_hat)%2 == 0):
                status = 0
                print("Parity check is satisfied, we can HALT after %d updates!\n" %(i+1))
                break
            else:
                print("Parity check is NOT satisfied after %d updates"%(i+1))
        if i == max_iter-1:
            print("Parity check cannot be satisfied after %d updates. Failed!\n" % max_iter)
        return x_hat, status
        
results = []
for i in range(31):
    ASC = 0
    for j in range(8):
        ASC += x_hat_1[8*i+j]*2**(7-j)
    results.append(str(chr(int(ASC))))
    print(x_hat_1[8*i:8*(i+1)],chr(int(ASC)))
''.join(results)
