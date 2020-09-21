#识别参数
Args <- commandArgs()
if(length(Args)==11){
  work_dir <- print(Args[6])
  tad_dir_file <- print(Args[8])
  mcool_dir_file <- print(Args[7])
  compartment_dir_file <- print(Args[9])
  out_dir <- print(Args[10])
  chr_num <- print(Args[11])
}
if(length(Args)==12){
  work_dir <- print(Args[6])
  tad_dir_file <- print(Args[8])
  mcool_dir_file <- print(Args[7])
  compartment_dir_file <- print(Args[9])
  out_dir <- print(Args[10])
  chr_num <- print(Args[11])
  com_name <- print(Args[12])
}
if(length(Args)==13){
  work_dir <- print(Args[6])
  tad_dir_file <- print(Args[8])
  mcool_dir_file <- print(Args[7])
  compartment_dir_file <- print(Args[9])
  out_dir <- print(Args[10])
  chr_num <- print(Args[11])
  com_name <- print(Args[12])
  tad_name <- print(Args[13])
}

#加载包
library("tidyr")
library("visNetwork")

#记录程序开始的时间
start_time <- Sys.time()
#设置工作路径
setwd(work_dir)
#读取对应染色体的.mcool文件
chr_mcool <- read.table(mcool_dir_file,header=T)
#读入compartment文件
all_compartment <- read.csv(compartment_dir_file,header=T,sep = "\t")
#提取对应染色体的compartment信息
chr_compartment <- all_compartment[which(all_compartment$chrom==chr_num),]
chr_compartment$com_name <- paste("com_",c(1:nrow(chr_compartment)),sep="")
#读取TAD位置信息
all_tad <- read.table(tad_dir_file,header=T)
#提取不同染色体的tad信息
chr_tad <- all_tad[which(all_tad$chrom == chr_num),]
chr_tad[,2:4] <- as.data.frame(lapply(chr_tad[,2:4],as.numeric))
chr_tad$tad_name <- paste("tad_",c(1:nrow(chr_tad)),sep="")
#去除端粒的交互信息
min <- which(chr_mcool$start1==chr_compartment$start[1])[1]
max <- tail(which(chr_mcool$end1==chr_compartment$end[nrow(chr_compartment)]),n=1)
chr_mcool <- chr_mcool[min:max,]
chr_mcool <- chr_mcool[order(chr_mcool$start2),]
row.names(chr_mcool) <- c(1:nrow(chr_mcool))
min <- which(chr_mcool$start2==chr_compartment$start[1])[1]
max <- tail(which(chr_mcool$end2==chr_compartment$end[nrow(chr_compartment)]),n=1)
chr_mcool <- chr_mcool[min:max,]
row.names(chr_mcool) <- c(1:nrow(chr_mcool))
chr_mcool <- chr_mcool[order(chr_mcool$start1),]
#compartment的交互网络
if(length(Args)==11){
#if(1==1){
  #统计compartment间的交互信息
  chrom1_com <- c()#用于存储mcool文件中每行chrom1所在compartment
  chrom2_com <- c()#用于存储mcool文件中每行chrom2所在compartment
  m <- 1
  #通过循环，将mcool中每条交联信息位置对应到compartment位置
  for(i in 1:nrow(chr_mcool)){
    tad1 <- which(chr_mcool$start1[i] >= chr_compartment$start & chr_mcool$end1[i] <= chr_compartment$end)
    tad2 <- which(chr_mcool$start2[i] >= chr_compartment$start & chr_mcool$end2[i] <= chr_compartment$end)
    chrom1_com[m] <- paste('com_',tad1,sep="")
    chrom2_com[m] <- paste('com_',tad2,sep="")
    m <- m+1
  }
  #生成links文件
  compartment_links <- data.frame(chrom1_com,chrom2_com,chr_mcool$count)
  compartment_links <- tidyr::unite(compartment_links,"com",chrom1_com,chrom2_com,sep="-")
  compartment_links <- aggregate(compartment_links[,"chr_mcool.count"],by=list(compartment_links[,"com"]),mean)
  compartment_links <- separate(compartment_links,Group.1,c("source","target"),sep="-",remove = T)
  com_count <- length(unique(chrom1_com))
  colnames(compartment_links) <- c("source","target","value")
  #过滤掉link数小于平均值的线
  compartment_links <- compartment_links[which(compartment_links$value>mean(compartment_links$value)),]
  #给compartment分组
  groups <- chr_compartment$com
  groups[which(chr_compartment$E1>0)] <- "compartment_A"
  groups[which(chr_compartment$E1<0)] <- "compartment_B"
  #计算nodes数
  #点的大小以边的count值的加和来计算
  value <- com_count
  for(i in 1:length(com_count)){
    value[i] <- sum(compartment_links[which(com_count[i]==compartment_links$source),]$value)
  }
  #生成nodes文件
  compartment_nodes <-data.frame(id=chr_compartment$com,
                                 label=chr_compartment$com,
                                 group=groups,
                                 title=paste(chr_compartment$start,chr_compartment$end,sep="_"),
                                 value=value)
  #生成edges文件
  compartment_edges <- data.frame(from=compartment_links$source,to=compartment_links$target)
  #作图并存储
  compartment_title <- paste("network_node_compartment_",chr_num,sep = "")
  compartment_network <- visNetwork(compartment_nodes,compartment_edges,main=compartment_title)%>%
    visOptions(selectedBy = "group")%>%
    visLegend(useGroups = TRUE,width = 0.3,position = "right")%>%
    visOptions(highlightNearest = TRUE)%>%
    visInteraction(dragNodes = TRUE,## 可移动节点
                   dragView = TRUE, # 移动图
                   zoomView = TRUE,  # 缩放
                   navigationButtons = TRUE # 下方添加按钮
    )%>%
    visLayout(randomSeed = 4)
  compartment_out_dir_file <- paste(out_dir,compartment_title,sep="")
  visSave(compartment_network,file = paste(compartment_out_dir_file,".html"))
}

#选取指定compartment中的所有tad进行绘图
if(length(Args)==12 || length(Args)==13){
  #将tad划分至compartment中
  m <- 1
  chr_tad$tad_com <- c(1:nrow(chr_tad))
  for(i in 1:nrow(chr_tad)){
    if(chr_tad$end[i] < chr_compartment$start[m]){
      next
    }
    else if(chr_tad$start[i] >= chr_compartment$start[m]){
      while(chr_tad$start[i] >= chr_compartment$end[m]){
        if(m<nrow(chr_compartment)){
          m <- m+1
        }
        else{
          break
        }
      }
      if(chr_tad$end[i] <= chr_compartment$end[m]){
        chr_tad$tad_com[i] <- chr_compartment$com[m]
      }
      else if(chr_tad$end[i] > chr_compartment$end[m] & m<nrow(chr_compartment)){
        chr_tad$tad_com[i] <- "2_com"
        if(m < nrow(chr_compartment))
        {
          m <- m+1
        }
        else{
          m <- m
        }
      }
    }
  }
  
  chr_tad_com <- chr_tad[which(chr_tad$tad_com==com_name),]
  if(is.na(chr_tad_com$chrom[1])){ #如果compartment中不含tad，打印以下语句
    print("No tad totally belong to this compartment!")
  }
  else{
    #提取tad内的compartment
    min_1 <- min(which(chr_mcool$start1 == chr_tad_com$start[1]))
    max_1 <- max(which(chr_mcool$end1 == chr_tad_com$end[nrow(chr_tad_com)]))
    chr_mcool_tad <- chr_mcool[min_1:max_1,]
    chr_mcool_tad <- chr_mcool_tad[order(chr_mcool_tad$start2),]
    min_2 <- min(which(chr_mcool_tad$start2 == chr_tad_com$start[1]))
    max_2 <- max(which(chr_mcool_tad$end2 == chr_tad_com$end[nrow(chr_tad_com)]))
    chr_mcool_tad <- chr_mcool_tad[min_2:max_2,]
    #统计tad间的交互关系
    chrom1_tad <- c()
    chrom2_tad <- c()
    m <- 1
    for(i in 1:nrow(chr_mcool_tad)){
      tad1 <- which(chr_mcool_tad$start1[i] >= chr_tad_com$start & chr_mcool_tad$end1[i] <= chr_tad_com$end)
      tad2 <- which(chr_mcool_tad$start2[i] >= chr_tad_com$start & chr_mcool_tad$end2[i] <= chr_tad_com$end)
      chrom1_tad[m] <- chr_tad_com$tad_name[tad1]
      chrom2_tad[m] <- chr_tad_com$tad_name[tad2]
      m <- m+1
    }
    #生成links
    tad_links <- data.frame(chrom1_tad,chrom2_tad,chr_mcool_tad$count)
    tad_links <- tidyr::unite(tad_links,"tad",chrom1_tad,chrom2_tad,sep="-")
    tad_links <- aggregate(tad_links[,"chr_mcool_tad.count"],by=list(tad_links[,"tad"]),mean)
    tad_links <- separate(tad_links,Group.1,c("source","target"),sep="-",remove = T)
    tad_count <- length(unique(chrom1_tad))
    colnames(tad_links) <- c("source","target","value")
    #点的大小以边的count值的加和来计算（效果不明显）
    value <- chr_tad_com$tad_name
    for(i in 1:length(chr_tad_com$tad_name)){
      value[i] <- sum(tad_links[which(chr_tad_com$tad_name[i]==tad_links$source),]$value)
    }
    #生成nodes文件
    tad_nodes <-data.frame(id=chr_tad_com$tad_name,
                           label=chr_tad_com$tad_name,
                           title=paste(chr_tad_com$start,chr_tad_com$end,sep="_"),
                           value=value)
    #生edge文件
    tad_edges <- data.frame(from=tad_links$source,to=tad_links$target)
    #作图并存储
    tad_title <- paste("network_node_tad_",chr_num,sep="")
    tad_title <- paste(tad_title,"_",sep="")
    tad_title <- paste(tad_title,com_name,sep="")
    tad_network <- visNetwork(tad_nodes,tad_edges,main=tad_title)%>%
      visInteraction(dragNodes = TRUE,## 可移动节点
                     dragView = TRUE, # 移动图
                     zoomView = TRUE,  # 缩放
                     navigationButtons = TRUE # 下方添加按钮
      )%>%
      visIgraphLayout(layout = "layout_in_circle")
    tad_out_dir_file <- paste(out_dir,tad_title,sep="")
    visSave(tad_network,file = paste(tad_out_dir_file,".html",sep=""))
    
    if(length(Args)==13){
      min_1 <- min(which(chr_mcool_tad$start1==chr_tad_com$start[which(chr_tad_com$tad_name==tad_name)]))
      max_1 <- max(which(chr_mcool_tad$end1==chr_tad_com$end[which(chr_tad_com$tad_name==tad_name)]))
      chr_mcool_tad_2 <-chr_mcool_tad[min_1:max_1,]
      chr_mcool_tad_2 <- chr_mcool_tad_2[order(chr_mcool_tad_2$start2),]
      min_2 <- min(which(chr_mcool_tad$start2==chr_tad_com$start[which(chr_tad_com$tad_name==tad_name)]))
      max_2 <- max(which(chr_mcool_tad$end2==chr_tad_com$end[which(chr_tad_com$tad_name==tad_name)]))
      chr_mcool_tad_2 <- chr_mcool_tad_2[min_2:max_2,]
      
      #统计fragement之间的交互信息
      fragment_links <- tidyr::unite(chr_mcool_tad_2,"fragment",c("start1","end1"),sep="_")
      fragment_links <- tidyr::unite(fragment_links,"fragment_2",c("start2","end2"),sep="_")
      fragment <- unique(fragment_links$fragment)
      fragment_name <- paste("fragment",c(1:length(fragment)))
      for(i in 1:length(fragment)){
        fragment_links$fragment_name[match(fragment[i],fragment_links$fragment)] <- fragment_name[i]
        fragment_links$fragment_2_name[match(fragment[i],fragment_links$fragment_2)] <- fragment_name[i]
      }
      #点的大小以边的count值的加和来计算
      fragment_value <- fragment
      for(i in 1:length(fragment)){
        fragment_value[i] <- sum(fragment_links[which(fragment_name[i]==fragment_links$fragment_name),]$count)
      }
      #生成fragment图的nodes和edges
      fragment_nodes <- data.frame(id=fragment_name,
                                   label=fragment_name,
                                   title=fragment,
                                   vlaue=fragment_value)
      fragment_edges <- data.frame(from=fragment_links$fragment_name,
                                   to=fragment_links$fragment_2_name,
                                   value=(fragment_links$count)/100)
      #画fragement网络图并储存
      fragment_title <- paste(tad_title,"_",sep="")
      fragment_title <- paste(fragment_title,tad_name,sep="")
      network <- visNetwork(fragment_nodes,fragment_edges,main=fragment_title)%>%
        visInteraction(dragNodes = TRUE,## 可移动节点
                       dragView = TRUE, # 移动图
                       zoomView = TRUE,  # 缩放
                       navigationButtons = TRUE # 下方添加按钮
        )%>%
        visLayout(randomSeed = 4)
      fragment_out_dir_file <- paste(out_dir,fragment_title,sep="")
      visSave(network,file = paste(fragment_out_dir_file,".html",sep=""))
    }
  }
}

#记录程序结束的时间
end_time <- Sys.time()
run_time <- end_time - start_time
print(paste("running time: ",run_time))