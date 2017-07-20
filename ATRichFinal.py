
# coding: utf-8

# In[1]:


import numpy as np
dict_ = {'A':0, 'T':1, 'G':2,'C':3}


# In[2]:


p = open("/home/vatsala/Desktop/DeepBind/nbt3300-code/data/dream5/chipseq/TF_23_CHIP_51_full_genomic.seq")

x_train=[]
y_train=[]
temp = []


# In[3]:


for line in p:
    temp.append(line)

#print temp


# In[4]:


arr=np.array(temp)
print arr[0][51]


# In[5]:


ans = np.zeros((51,4))
for ox in range(52):
    for ix in range(1000):
        if arr[ix][ox]=='A':
            ans[ox][0]+=1

        if arr[ix][ox]=='T':
            ans[ox][1]+=1

        if arr[ix][ox]=='G':
            ans[ox][2]+=1

        if arr[ix][ox]=='C':
            ans[ox][3]+=1
    

#print ans

prob = (ans/1000.0)
#print prob


maxArr = np.zeros((50,1))




for ix in range(50):
    store=max(prob[ix][0],prob[ix][1],prob[ix][2],prob[ix][3])
    maxArr[ix][0]=store





print maxArr
#np.save('Consensus_Array',ans)


# In[6]:


import seaborn as sns
get_ipython().magic(u'matplotlib inline')
#uniform_data = np.random.rand(10, 12)
ax = sns.heatmap(maxArr, vmin=0.20, vmax=0.40)


# In[7]:


ax=sns.heatmap(prob, xticklabels='ATGC')


# In[8]:


flank1 = np.zeros((50,4))
for ox in range(10):
    for ix in range(1000):
        if arr[ix][ox]=='A':
            flank1[ox][0]+=1

        if arr[ix][ox]=='T':
            flank1[ox][1]+=1

        if arr[ix][ox]=='G':
            flank1[ox][2]+=1

        if arr[ix][ox]=='C':
            flank1[ox][3]+=1


for ox in range(40,50):
    for ix in range(1000):
        if arr[ix][ox]=='A':
            flank1[ox][0]+=1

        if arr[ix][ox]=='T':
            flank1[ox][1]+=1

        if arr[ix][ox]=='G':
            flank1[ox][2]+=1

        if arr[ix][ox]=='C':
            flank1[ox][3]+=1
    

print flank1
np.save('flank1_Array',flank1)


# In[9]:


ax=sns.heatmap(flank1, xticklabels='ATGC')


# In[10]:


maxArrWithBaseInfo = np.zeros((50,2))
#print maxArrWithBaseInfo
Index=np.argmax(prob,axis=1)
print Index
for ix in range(50):
    store=max(prob[ix][0],prob[ix][1],prob[ix][2],prob[ix][3])
    
    
    maxArrWithBaseInfo[ix][0]=store
    maxArrWithBaseInfo[ix][1]=Index[ix]
    
#print maxArrWithBaseInfo

seq=""
for ix in range(51):
    if Index[ix]==0:
        seq+='A'
    if Index[ix]==1:
        seq+='T'
    if Index[ix]==2:
        seq+='G'
    if Index[ix]==3:
        seq+='C'
print seq


# In[11]:


minorStartIndex=np.zeros([51])
numberOfATrichSites=0
j=0
for i in range(50):
    if(seq[i]=='A'and seq[i+1]=='T'):
        numberOfATrichSites+=1;
        minorStartIndex[j]=i;
        j+=1;
print ("Number of AT rich sites are = "+str(numberOfATrichSites))
print ("Indices of AT rich sequences are:")

for i in range(numberOfATrichSites):
                           print minorStartIndex[i], minorStartIndex[i]+1

        

