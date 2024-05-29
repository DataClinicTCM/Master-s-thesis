#######################################
### Project: Network proximity calculation.
### Author: 参考Deisy Gysi的代码，使用了Netsci包
### date: March 5st 2024
#######################################

rm(list = ls())
library(tidyr)
library(magrittr)
library(data.table)
library(igraph)
library(dplyr)
library(NetSci)

# 由于源代码的疾病是COVID-19,后面都是用这个称呼其他疾病(懒得改名称)

PPI = fread("DatasetS2.csv") # Read PPI
COVID = fread("disease_targets.csv") # Read COVID Genes
PPI %<>% filter(!is.na(proteinA_entrezid) & !is.na(proteinB_entrezid))
gPPI = graph_from_data_frame(PPI, directed = F) %>% 
  igraph::simplify() # Create the PPI network and remove self-loops

# ----------------------------------------------------
# V(g)/V(g)$label/V(g)$name不一样!!!!!!!!!
# V(g)返回索引，是igraph对象
# V(g)$label不指定的话通常为空，返回NULL(与NA不同)
# V(g)$name返回节点名称，是character
# 由于%in%运行时会进行隐式的类型转换，所以要提前检查
# ----------------------------------------------------
# 判断疾病靶点是否在PPI中,first judgment data type
class(COVID$ENTREZID)
class(V(gPPI)$name)
class(V(gPPI))

as.character(COVID$ENTREZID) %in% V(gPPI)$name
sum(as.character(COVID$ENTREZID) %in% V(gPPI)$name)

# 选择在PPI中的靶点
COVID$ENTREZID[as.character(COVID$ENTREZID) %in% V(gPPI)$name]
class(COVID$ENTREZID[as.character(COVID$ENTREZID) %in% V(gPPI)$name])

# Target_COVID = COVID %>% 
#   pull(ENTREZID) %>% 
#   as.character()

Target_COVID <- COVID$ENTREZID[as.character(COVID$ENTREZID) %in% V(gPPI)$name]
class(Target_COVID)

gCOVID = gPPI %>% 
  induced_subgraph(., vids = as.character(Target_COVID)) # Select the subgraph from COVID

# Plot it
V(gCOVID)$label = NA
gCOVID %<>% simplify() #去除自连接和多条边
V(gCOVID)$size = degree(gCOVID) #按照degree计算节点大小
V(gCOVID)$size = degree(gCOVID) %>% CoDiNA::normalize() #节点大小归一化
V(gCOVID)$size  = (V(gCOVID)$size + 0.1)*5 #调整大小
gCOVID %<>% delete_vertices(., degree(.) == 0) #删除degree=0的节点
coord = layout_with_fr(gCOVID) #力导向布局
V(gCOVID)$label = ifelse(V(gCOVID)$size > 1.5, V(gCOVID)$name, NA) #节点大于2.5才有名字

# 将名字由ENTREZID改为SYMBOL(这段可注释掉)
# V(gCOVID)$label
aa <- V(gCOVID)$label[which(!is.na(V(gCOVID)$label))]
dput(V(gCOVID)$label[which(!is.na(V(gCOVID)$label))])
eg1 = clusterProfiler::bitr(aa, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
# V(gCOVID)$label[which(!is.na(V(gCOVID)$label))] <- eg1$SYMBOL

# =============
# 把eg1作为后面GO KEGG的gene输入
# ==============



E(gCOVID)$curved =0.1
V(gCOVID)$color = V(gCOVID)$frame.color = "red"
V(gCOVID)$label.color = "black"
V(gCOVID)$label.cex = V(gCOVID)$size/5
E(gCOVID)$color = 'salmon1'
plot(extract_LCC(gCOVID)) #plot LCC

# Calculate the LCC significance, degree preserving
set.seed(666)
LCC = LCC_Significance(N = 1000, 
                 Targets = as.character(Target_COVID), 
                 G = gPPI,
                 min_per_bin = 30)

# Plot the distribution
Histogram_LCC(LCC)

# ----------------------------------------------------------------------------
# use NetSci package to calculate proximity
# ----------------------------------------------------------------------------
# 最好设置随机种子
set.seed(666)
# Extract LCC
PPI_LCC = extract_LCC(gPPI)

# 计算3个药物之间的分离度Sab---------------------------------
# read data
SDH <- read.csv(file = "SDH_targets.csv")
SY <- read.csv(file = "SY_targets.csv")
SZY <- read.csv(file = "SZY_targets.csv")

# 判断药物靶点是否在PPI_LCC中
sum(as.character(SDH$ENTREZID) %in% V(PPI_LCC)$name)
sum(as.character(SY$ENTREZID) %in% V(PPI_LCC)$name)
sum(as.character(SZY$ENTREZID) %in% V(PPI_LCC)$name)


# 选择在PPI中的药物靶点
SDH$ENTREZID[as.character(SDH$ENTREZID) %in% V(PPI_LCC)$name]
class(SDH$ENTREZID[as.character(SDH$ENTREZID) %in% V(PPI_LCC)$name])
Target_SDH <- SDH$ENTREZID[as.character(SDH$ENTREZID) %in% V(PPI_LCC)$name]

SY$ENTREZID[as.character(SY$ENTREZID) %in% V(PPI_LCC)$name]
class(SY$ENTREZID[as.character(SY$ENTREZID) %in% V(PPI_LCC)$name])
Target_SY <- SY$ENTREZID[as.character(SY$ENTREZID) %in% V(PPI_LCC)$name]

SZY$ENTREZID[as.character(SZY$ENTREZID) %in% V(PPI_LCC)$name]
class(SZY$ENTREZID[as.character(SZY$ENTREZID) %in% V(PPI_LCC)$name])
Target_SZY <- SZY$ENTREZID[as.character(SZY$ENTREZID) %in% V(PPI_LCC)$name]


# 构建输入数据
D1 = data.frame(gene = Target_SDH, drug = "SDH")
D2 = data.frame(gene = Target_SY, drug = "SY")
D3 = data.frame(gene = Target_SZY, drug = "SZY")
Drugs = rbind(D1, D2, D3)
Drugs %<>% dplyr::select(drug, gene)

set.seed(666)
# calculate
separation_Significance(G = PPI_LCC,
                        ST = Drugs,
                        correct_by_target = FALSE,
                        Threads = 2)

# 计算疾病与药物的邻近度-----------------------------------------------------
########### 注意!!!!! set,ST数据名称必须要与ID列一致 #############
SDH <- D1
SY <- D2
SZY <- D3 #注意SDH SY SZY不是原始数据了

colnames(SDH) <- c("Target", "ID")
colnames(SY) <- c("Target", "ID")
colnames(SZY) <- c("Target", "ID")

SDH %<>% dplyr::select(ID, Target)
SY %<>% dplyr::select(ID, Target)
SZY %<>% dplyr::select(ID, Target)

# 顶点格式要转换为character
SDH$Target <- as.character(SDH$Target)
SY$Target <- as.character(SY$Target)
SZY$Target <- as.character(SZY$Target)

set.seed(666)
# 计算
avr_proximity_multiple_target_sets(
  set = c("SDH", "SY", "SZY"),
  G = PPI_LCC,
  ST = rbind(SDH, SY, SZY),
  source = as.character(Target_COVID),
  N = 1000,
  bins = 100,
  min_per_bin = 20
) ###################### need long time,maybe 40 min



# ----------------------------------------------- output to Gephi,cytoscape
write.csv(Target_COVID, file = "diseaseID.csv")

write.csv(Target_SDH, file = "SDH-ID.csv")
write.csv(Target_SY, file = "SY-ID.csv")
write.csv(Target_SZY, file = "SZY-ID.csv")

# ----------------------------------------------------后续完善该部分--------
# 当获取到在PPI中的药物靶点之后，就可以直接把这三类靶点加到gCOVID网络中，就能够画出来
# 然后以不同颜色突出显示这三类靶点

gOS_SSS = gPPI %>% 
  induced_subgraph(., vids = as.character(unique(c(
                               Target_COVID,
                               Target_SDH,
                               Target_SY,
                               Target_SZY))))


plot(gOS_SSS) #ugly

# plot it beautiful
V(gOS_SSS)$label = NA
gOS_SSS %<>% simplify() #去除自连接和多条边
V(gOS_SSS)$size = degree(gOS_SSS) #按照degree计算节点大小
V(gOS_SSS)$size = degree(gOS_SSS) %>% CoDiNA::normalize() #节点大小归一化
V(gOS_SSS)$size  = (V(gOS_SSS)$size + 0.1)*5 #调整大小
gOS_SSS %<>% delete_vertices(., degree(.) == 0) #删除degree=0的节点
# coord1 = layout_in_circle(gOS_SSS) #貌似没啥用
V(gOS_SSS)$label = ifelse(V(gOS_SSS)$size > 2.5, V(gOS_SSS)$name, NA) #节点大于2.5才有名字

# 将名字由ENTREZID改为SYMBOL(这段可注释掉)
# V(gOS_SSS)$label
# aa <- V(gOS_SSS)$label[which(!is.na(V(gOS_SSS)$label))]
# dput(V(gOS_SSS)$label[which(!is.na(V(gOS_SSS)$label))])
# eg2 = clusterProfiler::bitr(aa, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
# V(gOS_SSS)$label[which(!is.na(V(gOS_SSS)$label))] <- eg2$SYMBOL

E(gOS_SSS)$curved =0.1
V(gOS_SSS)$color = V(gOS_SSS)$frame.color = "red"
V(gOS_SSS)$label.color = "black"
V(gOS_SSS)$label.cex = V(gOS_SSS)$size/5
E(gOS_SSS)$color = 'salmon1'
plot(gOS_SSS, layout = layout_with_fr)

# 给不同靶点上色
# 画出来很丑，没办法实现不同靶点一个颜色并加圈，可以尝试画模式图


# ===============================================================================
# GO KEGG
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

library(fanyi)
set_translate_option(appid = "******************",
                     key = "********************",
                     source = "baidu")


targets <- readxl::read_xlsx(path = "cytoscape/network_compounds.xlsx", sheet = 3)
targets <- targets$target

#基因名称转换，返回的是数据框
gene = bitr(targets, 
            fromType="SYMBOL", 
            toType="ENTREZID", 
            OrgDb="org.Hs.eg.db")

# ===========
gene <- eg1 #对应于LCC大于1.5的点的富集分析
# ===========

# KEGG
EGG <- enrichKEGG(gene = gene$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)

#从KEGG富集分析结果中提取结果表格：
KEGG_result <- EGG@result

#先提取富集结果表前Top30：
KEGG_top30 <- KEGG_result[1:30,]

#指定绘图顺序（转换为因子）：
KEGG_top30$pathway <- factor(KEGG_top30$Description,levels = rev(KEGG_top30$Description))

#Top30显著性气泡图：
#将pathway按照p值排列：
p1 <- ggplot(data = KEGG_top30,
       aes(x = Count,
           y = pathway))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+ # 气泡大小及颜色设置
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched KEGG Pathways",
       size = "Count")
# --------------------------------
# Visualize enriched KEGG pathways
head(EGG) #first pathway is wnt
browseKEGG(EGG, 'hsa04310')

data(geneList, package="DOSE") #DOSE genelist里面带值，但是与本研究好像无关，所以不要用颜色映射的
# 如果有自己的genelist可以做
library("pathview")
hsa04310 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04310",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
# --------------------------------
# pmcplot(KEGG_top30$Description[1], 2020:2024, proportion=FALSE)
# --------------------------------

# GO
EGO <- enrichGO(gene = gene$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

#从GO富集分析结果中提取结果表格：
GO_result <- EGO@result

#先提取富集结果表各分类前10：
# GO_top30 <- GO_result %>%
#   group_by(ONTOLOGY) %>%
#   filter(row_number() <= 10) #注意filter用法
# 或
GO_top30 <- GO_result %>% 
  group_by(ONTOLOGY) %>% 
  dplyr::slice_head(n = 10) #前10行


#指定绘图顺序（转换为因子）：
GO_top30$Description <- factor(GO_top30$Description,levels = rev(GO_top30$Description))

#Top30显著性气泡图：
#将Description按照p值排列：
p2 <- ggplot(data = GO_top30,
       aes(x = Count,
           y = Description))+
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+ # 气泡大小及颜色设置
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  facet_grid(ONTOLOGY~., scale = "free") +
  labs(x = "Gene Number",
       y = "",
       title = "Dotplot of Enriched GO",
       size = "Count")

# 用Y叔的包汉化(记得赋值，不然下次画图可能还用额度)
# p3 <- translate_ggplot(p1, axis="y")
# p4 <- translate_ggplot(p2, axis="y")
# aplot::plot_list(KEGG = p1, GO = p2, ncol=2)

# ---------------treeplot----------------------
EGOx <- enrichplot::pairwise_termsim(EGO)
enrichplot::treeplot(EGOx)

# ==========================Comparing multiple gene lists=======================

sdh <- readxl::read_xlsx(path = "结果整合.xlsx", sheet = 2)
szy <- readxl::read_xlsx(path = "结果整合.xlsx", sheet = 3)
sy <- readxl::read_xlsx(path = "结果整合.xlsx", sheet = 4)

sdh_xcs <- sdh[sdh$compound == "香草酸",]$gene
sdh_dhkgy <- sdh[sdh$compound == "地黄苦苷元",]$gene
sdh_jnpgs <- sdh[sdh$compound == "京尼平苷酸",]$gene
sdh_zzxg <- sdh[sdh$compound == "栀子新苷",]$gene

szy_mszs <- szy[szy$compound == "没食子酸",]$gene
szy_mng <- szy[szy$compound == "莫诺苷",]$gene
szy_yecs <- szy[szy$compound == "原儿茶酸",]$gene
szy_zycg <- szy[szy$compound == "獐牙菜苷",]$gene

sy_gec <- sy$gene

gc <- list()
dput(ls())
compound_target <- c("sdh_dhkgy", "sdh_jnpgs", "sdh_xcs", "sdh_zzxg",
                     "sy_gec",
                     "szy_mng", "szy_mszs", "szy_yecs", "szy_zycg")

# 循环转换ID
for (i in compound_target) {
  gc[[i]] = bitr(eval(parse(text = i)), 
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")$ENTREZID
}

names(gc) <- c("sdh_地黄苦苷元", "sdh_京尼平苷酸", "sdh_香草酸", "sdh_栀子新苷",
               "sy_庚二醇",
               "szy_莫诺苷", "szy_没食子酸", "szy_原儿茶酸", "szy_獐牙菜苷")

# GO
ck <- compareCluster(geneCluster = gc, #a list of entrez gene id
                     fun = enrichGO, #enrich type
                     OrgDb='org.Hs.eg.db')

# ck <- pairwise_termsim(ck)
# ck_sim <- simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)

pgo <- dotplot(ck, font.size = 10, label_format = 60) + labs(x = NULL)+
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
print(pgo)



pgo_zh <- translate_ggplot(pgo, axis="y")
print(pgo_zh)

# ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# cnetplot(ck)

# KEGG
ck1 <- compareCluster(geneCluster = gc, #a list of entrez gene id
                     fun = enrichKEGG)

pgo1 <- dotplot(ck1, font.size = 10, label_format = 60) + labs(x = NULL)+
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
print(pgo1)


pgo_zh1 <- translate_ggplot(pgo1, axis="y")
print(pgo_zh1)

# --------------------pdf show chinese-------------------------
library(showtext)
font_add("songti","simsun.ttc")
font_add("TNM","times.ttf")
showtext_auto(enable=TRUE)

pgo <- dotplot(ck, font.size = 10, label_format = 60) + labs(x = NULL)+
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        # axis.text = element_text('TNM'), 
        axis.title = element_text("songti"))
ggsave(filename = 'pgo.pdf', plot = pgo, width = 8, height = 8)



ck1 <- compareCluster(geneCluster = gc, #a list of entrez gene id
                      fun = enrichKEGG)

pkegg <- dotplot(ck1, font.size = 10, label_format = 60) + labs(x = NULL)+
  scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        # axis.text = element_text('TNM'), 
        axis.title = element_text("songti"))
ggsave(filename = 'pkegg.pdf', plot = pkegg, width = 8, height = 8)
