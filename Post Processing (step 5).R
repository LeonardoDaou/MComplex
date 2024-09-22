Post_Process=function(){
  library(foreach)
  library(doParallel)
  
  ranking=function(cluster,proteins,g){
    returny=as.data.frame(matrix(0,ncol=2,nrow=length(proteins)))
    returny[,1]=proteins
    colnames(returny)=c("Proteins","Degree")
    for(i in 1:length(proteins)){
      returny[i,2]=sum(g[cluster,i])
    }
    returny=returny[order(returny$Degree,decreasing = TRUE),]$Proteins
  }
  data=function(path,sep=sep){
    collins_ppi=read.table(path,sep=sep)
    collins=unique(c(as.vector(collins_ppi[,2]),as.vector(collins_ppi[,1])))
    s=matrix(nrow = length(collins),ncol=length(collins))
    colnames(s)=collins
    rownames(s)=collins
    i=1
    while(i<=length(collins_ppi[,2])){
      s[collins_ppi[i,1],collins_ppi[i,2]]=1#collins_ppi[i,3]#or 1
      s[collins_ppi[i,2],collins_ppi[i,1]]=1#collins_ppi[i,3]#or 1
      i=i+1
    }
    s[which(is.na(s))]=0
    return(s)
  }
  cohesiveness=function(cluster,g,p){
    ins=sum(g[cluster,cluster])/2
    locs=vector()
    for(i in 1:length(cluster)){
      locs=c(locs,which(colnames(g)==cluster[i]))
    }
    outs=sum(g[cluster,-locs])
    coh=ins/(outs+ins+p)
    return(coh)
  }
  
  degs_rem=function(clusters,g){
    returny=vector(length=1,mode="list")
    rems=vector(length=1,mode="list")
    
    mat=as.data.frame(matrix(ncol=3))
    colnames(mat)=c("complex","ext_deg","length")
    current=clusters
    remaining=unlist(current[which(lengths(current)==1)])
    if(is.null(remaining)==FALSE){
      current=current[-which(lengths(current)==1)]
    }
    for(j in 1:length(current)){
      current_cl=current[[j]]
      len=length(current_cl)
      outs=NA
      if(len>1){
        locs=vector()
        inter=intersect(remaining,current_cl)
        copied=setdiff(remaining,inter)
        outs=length(unique(which(g[current_cl,copied]>0)))
      }
      m=as.data.frame(matrix(nrow=1,ncol=3))
      m[1,1]=list(current[j])
      m[1,2]=outs
      m[1,3]=len
      colnames(m)=c("complex","ext_deg","length")
      mat=rbind(mat,m)
    }
    mat=mat[-1,]
    returny[[1]]=mat
    rems[[1]]=remaining
    
    rets=vector(length=2,mode="list")
    rets[[1]]=returny
    rets[[2]]=rems
    names(rets)=c("exts","remainings")
    return(rets)
  }
  
  Cl_refine_one=function(clusters,runs=1,data,times=1){
    returny=vector(length=1,mode="list")
    for(k in 1:times){
      #print(paste("Timestep-----------------------------------------",k,sep=""))
      returny_sub=vector(length=length(clusters),mode="list")
      q=clusters
      for(j in 1:runs){
        print("Processing the degrees")
        deggy=degs_rem(q,data)
        exts=deggy$exts
        remainingsy=deggy$remainings
        new_q=vector(mode="list",length=1)
        a=1
        for(i in 1:a){
          #print(paste("RUN: ",j,"/",runs," SET: ",i,"/",length(q),sep=""))
          print("Starting Merging Process")
          new_q[[i]]=vector(mode="list")
          current=clusters
          clusterings=vector(mode="list")
          if(length(remainingsy)!=0){
            remaining=remainingsy[[i]]
          }else{
            remaining=remainingsy
          }
          exts[[i]]=exts[[i]][order(exts[[i]]$ext_deg,decreasing = TRUE),]
          clus=exts[[i]]$complex
          for(l in 1:length(clus)){
            cl=clus[[l]]
            copy_cl=cl
            for(ps in copy_cl){
              locs=vector()
              for(p in 1:length(cl)){
                locs=c(locs,which(colnames(data)==cl[p]))
              }
              outs=length(which(data[ps,-locs]>0))
              if(length(cl)>1){
                if(outs>1){
                  old_c=cohesiveness(cl,data,0)
                  new_cl=cl[-which(cl==ps)]
                  new_c=cohesiveness(new_cl,data,0)
                  if(new_c>old_c){
                    cl=new_cl
                    remaining=c(remaining,ps)
                  }
                }
              }
            }
            if(length(cl)>0){
              clusterings[[length(clusterings)+1]]=cl
            }
          }
          newer=vector(mode="list")
          newer=c(clusterings,as.list(remaining))
          deggies=degs_rem(newer,data)
          remaining=vector()
          if(length(deggies$remainings)>0){
            remaining=deggies$remainings[[1]]
          }
          exties=deggies$exts[[1]][order(deggies$exts[[1]]$ext_deg,decreasing = TRUE),]
          clus=exties$complex
          for(l in 1:length(clus)){
            cl=clus[[l]]
            if(length(remaining)>0 && !is.na(remaining[1])){
              inter=intersect(remaining,cl)
              compared=setdiff(remaining,inter)
              if(!is.null(dim(which(data[cl,compared]>0,arr.ind = TRUE)))){
                locs=unique(which(data[cl,compared]>0,arr.ind = TRUE)[,2])
                rems=colnames(data[cl,compared])[locs]
                old_coh=0
                new_coh=0
                if(length(rems)>0&&!is.na(rems[1])){
                  rems=ranking(cl,proteins = rems,g=data)
                  rs=1
                  while(rs<=length(rems)&&!is.na(rems[1])){
                    if(length(rems)==0){
                      break
                    }
                    old_coh=cohesiveness(cl,data,0)
                    new_coh=cohesiveness(c(cl,rems[rs]),data,0)
                    if(new_coh>old_coh){
                      cl=c(cl,rems[rs])
                      remaining=remaining[-which(remaining==rems[rs])[1]]
                      rems=rems[-rs]
                      rems=ranking(cl,proteins = rems,g=data)
                    }else{
                      rs=rs+1
                    }
                  }
                }
              }
            }
            new_q[[i]][[length(new_q[[i]])+1]]=cl
            
          }
          # go through different runs and contract and expand depending
          # if there is an increase or decrease of cohesion for the cluster
          # Create another one where density matters
          
          new_q[[i]]=c(new_q[[i]],as.list(remaining))
          print("End of Merging Process")
          save(new_q,file="returnies.RData")
          return(new_q)
        }
        q=new_q
      }
      #returny[[k]]=new_q
      
    }
  }
  merge_overlap_one=function(clusters,overper=0.8){
    current=clusters
    m=matrix(0,nrow=length(current),ncol=length(current))
    for(i in 1:length(current)){
      for(j in 1:length(current)){
        if(i!=j){
          over=overlap(current[[i]],current[[j]])
          if(over>overper){
            m[i,j]=over
          }
        }
      }
    }
    
    rownames(m)=1:length(current)
    colnames(m)=1:length(current)
    
    g=graph.adjacency(m,mode = "undirected",weighted = TRUE)
    connected=components(g)
    before=which(connected$csize>1)
    removing=c()
    temp=current
    print(length(before))
    print(before)
    for(i in before){
      a=which(connected$membership==i)
      cur=current[a]
      q=vector()
      for(l in cur){
        q=c(q,l)
      }
      q=unique(q)
      removing=c(removing,a)
      temp[[length(temp)+1]]=q
    }
    temp=temp[-removing]
    return(temp)
  }
  
  overlap=function(A,B){
    top=length(intersect(A,B))
    bottom=length(A)*length(B)
    return((top^2)/bottom)
  }
  
  dens=function(cluster,g){
    indices=vector(length=length(cluster))
    for(p in 1:length(cluster)){
      indices[p]=which(colnames(g)==cluster[p])
    }
    summy=sum(g[indices,indices])/2
    return(2*summy/(length(cluster)*(length(cluster)-1)))
  }
  
  clussy_each=function(locs,num){
    clusters=vector(mode="list",length=120)
    for(j in 1:6){
      print(paste("J--------------------",j))
      for(i in 1:num){
        if(j==1){
          clusters[[i]]=list()
        }
        print(paste("I--------",i))
        load(paste(locs,"/",j,"/All_Mapper_umap_p_",i,".RData",sep=""))
        a=gsub(perl=TRUE,pattern="V[0-9]+: ",replacement = "",e$Nodename)
        clus=unique(strsplit(x=a,split=", "))
        clusters[[i]]=unique(c(clusters[[i]],clus))
      }
    }
    return(clusters)
  }
  
  og_data=data(path = "./PPI Datasets/Gavin.txt",sep="\t")
  
  dir.create("Refined")
  clusters_all_pr_de_ec_label=clussy_each("./",num = 120)
  currently=vector(mode="list",length=120)
  totalCores=detectCores()
  clustery <- makeCluster(totalCores[1]-4)
  registerDoParallel(clustery)
  foreach(i=1:length(clusters_all_pr_de_ec_label))%dopar%{
    print(paste("Current ",i))
    library(igraph)
    processed_all_pr_de_ec_label=Cl_refine_one(clusters_all_pr_de_ec_label[[i]],data = og_data)[[1]]
    currently[[i]]=processed_all_pr_de_ec_label
    save(processed_all_pr_de_ec_label,file=paste("Refined/post_",i,".RData",sep=""))
  }
  stopCluster(clustery)
  
  currently=vector(mode="list",length=120)
  for(i in 1:120){
    load(paste("Refined/post_",i,".RData",sep=""))
    currently[[i]]=processed_all_pr_de_ec_label
  }
  dir.create("Post_Processed")
  values=c(0.5) #can add more if a wider is needed to test
  for(ths in values){
    dir.create(paste("Post_Processed/Overlap_",ths))
    clustery <- makeCluster(totalCores[1]-4)
    registerDoParallel(clustery)
    foreach(i=1:length(currently))%dopar%{
      library(igraph)
      processed_all_pr_de_ec_label_over=merge_overlap_one(currently[[i]],overper=ths)
      processed_all_pr_de_ec_label_no1=rem_1s(processed_all_pr_de_ec_label_over,data=s,th=0.2)
      processed_all_pr_de_ec_label_no1=unique(processed_all_pr_de_ec_label_no1$totality[[1]])
      for(j in 1:length(x)){
        write(paste(processed_all_pr_de_ec_label_no1[j][[1]],collapse = "\t"),file=paste("Post_Processed/Overlap_",ths,"/complexes_",j,".txt",sep=""),append=TRUE)
      }
    }
    stopCluster(clustery)
  }
}