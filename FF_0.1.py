####for fast filter
##Func. F_ff(or "importance" part): Fast-filtering function for segments. if the area decreases more than 20% by lossing child nodes (calculate by template), the segment will be defined as inactivate. This function is also used in checking procedure to ensure the anchored child nodes meet the minimum requirement.
##Func. ff-absence: Additional fast-filtering function for shadow nodes. It differs according to the type of shadow nodes.
####
##"node_list" is a dictionary that records every node in the netwrok.
##"centermap" is a (x,y,n)-sized list that records activated nodes. x: height of the input sample; y: width of the input sample; n: flexible length, records activated nodes in corresponding coordinate.

def importance_all(seg):
    #In a segment, for every combination of child nodes, if the lost area reaches 20%, the combination will be added to 'wontact'. The lost pure blocks will be ignored. In "wontact", every child node is IDed by its sequence in the template (key 'kernel' in the dictionary of a segment)
    
    if node_list[seg]['category']==1: #shadow nodes performing 'or' logic
        node_list[seg]['wontact']=[]
        return
    ll=len(node_list[seg]['kernel'])
    wontact=[set(range(ll))]
    base=[]
    for k in node_list[seg]['kernel']:
        base.append([k,1])
    #print(base)
    pure=[]
    od=list(range(ll))
    for i in range(ll):
        kname=node_list[seg]['kernel'][i][0]
        if 'pure' in node_list[kname].keys() and 'pure' not in node_list[seg].keys():
            od.remove(i)
            continue
        corr=copy.deepcopy(base)
        corr[i][1]=0
        flag=importancearea(seg,corr)
        if not flag:
            wontact.append(set([i]))
    
    comset=[]            
    for length in range(2,len(od)):
        data=[od]
        while len(data)<length:
            nextd=copy.deepcopy(data[-1])
            nextd.pop()
            data.append(nextd)
        combines=itertools.product(*data)#(range(ll),repeat=length)
        combines=list(combines)
        for item in combines:
            if len(set(item))<length:
                continue
            a=[]
            for i in item:
                kname=node_list[seg]['kernel'][i][0]
                if 'pure' in node_list[kname].keys() and 'pure' not in node_list[seg].keys():
                    continue
                a.append(i)
            if not a:
                continue
            #print('mark1',a)
            b=set(a)
            if b not in comset:
                comset.append(b)                                        
    #print('?!?!?!')
    for c in comset:
        corr=copy.deepcopy(base)
        for i in c:
            corr[i][1]=0
        flag=importancearea(seg,corr)
        
        if not flag:
            wontact.append(c)
    node_list[seg]['wontact']=copy.deepcopy(wontact)          
    #print('debug over')
    
def importancearea(seg,corres):
    #"seg" is the ID of a segment; "corres" records if a child node is considered "activated" (1 for activated and 0 for inactivated). The left area is calculated by the area covered by the activated child nodes (assume they are activated as expected). If the left area >=0.8 (which means the absent area<0.2), the combination is considered as "possible to activate the segment"; otherwise, the combination of inactivated child nodes will be recored in "wontact".
    print(seg)
    si,sj=node_list[seg]['center']
    area=[]
    count=0
    notpure=0
    #category=1 could use direct pure and notpure+1
    if node_list[seg]['category']==0:    
        for k in corres:
            if 'pure' not in node_list[k[0][0]].keys() or 'pure' in node_list[seg].keys():
                notpure+=1
                if not k[1]:
                    count+=1
                    continue
            else:
                continue
            if not k[1]:
                continue
            ci,cj=k[0][1],k[0][2]
            
            if node_list[k[0][0]]['level']==0:                
                for i,j in node_list[k[0][0]]['area']:
                    ind=[si+ci+i,sj+cj+j]
                    #ind=[k[1][1]+i,k[1][2]+j]
                    #ind=[si+i,sj+j]
                    if ind not in area:
                        area.append(ind)
            else:
                backarea=backseg_ini(k[0],seg)
                for ind in backarea:
                    #ind=[si+ci+i,sj+cj+j]
                    #ind=[i,j]
                    if ind not in area:
                        area.append(ind)
        if count>notpure*0.5:#int(notpure*(1-dthred)+0.5):#count>int(len(corres)*(1-dthred)+0.5):
            return False
        if len(area)<len(node_list[seg]['oriareanopure'])*0.8 and 'pure' not in node_list[seg].keys():#min(0.75,(notpure-1)/notpure):
            #(notpure-1)/notpure:#(maxsegblocknum-1)/maxsegblocknum: #*dthred*0.9375:#0.9375:#*0.95:#*0.9375:#*0.9375:#*0.875:#*dthred:
            return False
        if 'pure' in node_list[seg].keys() and len(area)<len(node_list[seg]['oriarea'])*0.8:
            return False
                        
    elif node_list[seg]['category']==1:#shadow C1
        print('?')
    
    return True        
    
def F_ff_absence(seg):
    #This function is part of the segment-calculating function. We list it seperatly to make the program more organized.
    k_syn=copy.deepcopy(node_list[syn]['kernel'])  
    wrongtolerance=len(k_syn)-int(len(k_syn)*dthred+0.5)#int(len(k_syn)*(1-dthred)+0.5)#int(len(k_syn)*0.1+0.5)
    k_kernel={}
    for k in k_syn:
        if k[0] not in k_kernel.keys():
            k_kernel[k[0]]=[]
        k_kernel[k[0]].append([k[1],k[2]])
        
    ks=[]
    for k in k_syn:#to count the kernel/seg needed to activate this seg
        ks.append(k[0])
    alen=len(ks)
    regks=dict(Counter(ks))#count the kernel/seg needed to activate this seg
    #for quickly decide if to calculate the seg or to return false
    copycmp=copy.deepcopy(centermap) 
    
    i,j=0,0#,max(1,int(wwi/4)),max(1,int(wwj/4)),max(1,int(wwi/4)),max(1,int(wwj/4))#i/x:vertical,j/y:horizenal
    forcalc_dic=[]
    acted_dic=[]
    acted_center=[]
    ##for every window:
        realarea=[]#calculate directly, no need to double-search the list#actual covered area
        ava={}
        for kk in k_kernel.keys():
            ava[kk]=0
        for di,dj in node_list[syn]['zoompadarea']:
            ri,rj=i+di,j+dj
            if ri>=0 and ri<inputshape[0] and rj>=0 and rj<inputshape[1]:
                realarea.append([ri,rj])
                #quick filter if this area contains enough syns
                #change in 2024.07
                for actk in copycmp[ri][rj]:
                    if actk in ava.keys():
                        ava[actk]+=1
        if not realarea:
            return False #area invalid
        irange=[0,0]
        jrange=[0,0]
        realarea.sort(key=lambda item: item[0])
        irange[0],irange[1]=realarea[0][0],realarea[-1][0]+1
        realarea.sort(key=lambda item: item[1])
        jrange[0],jrange[1]=realarea[0][1],realarea[-1][1]+1
        realarea.clear()
        for ri in range(irange[0],irange[1]):
            for rj in range(jrange[0],jrange[1]):
                realarea.append([ri,rj])
        ###editted in 2024.07
        ifcontinueflag=checkifcontinue(syn,ava,regks,k_kernel,wrongtolerance)
        if not ifcontinueflag:
            return False #skip current window            
        return True #calculate current window
        
def checkifcontinue(syn,ava,regks,k_kernel,wrongtolerance):
    #check if the missing nodes make it impossible to activate the segment
    losscount={}
    kk=copy.deepcopy(node_list[syn]['kernel'])
    loccount=0
    for k in k_kernel.keys():
        losscount[k]=[max(0,regks[k]-ava[k]),[]]
        if not losscount[k][0]:
            continue
        loccount+=losscount[k][0]
        for i in range(len(kk)):
            if kk[i][0]==k:
                losscount[k][1].append(i)
    if 'wontact' not in node_list[syn].keys():
        if loccount>wrongtolerance:
            return False
        return True
    data=[]
    for k in losscount.keys():
        for i in range(losscount[k][0]):
            data.append(losscount[k][1]) 
    
    u=itertools.product(*data)
    u=list(u)
    if () in u:
        u.remove(())
    combines=[]
    
    if u:
        for i in u:
            a=set(i)
            if len(a)<len(i):
                continue
            if a not in combines:
                combines.append(a)
    count=0
    for i in combines:
        if i not in node_list[syn]['wontact']:
            count+=1
            
    if not u:
        #print(11111)
        count=1
    if count>0:
        return True#further calcu
    else:
        return False
        

    
