library(vegan)
library(dplyr)

geodistancematrix<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/matrix/geo_distance10pop_matrix.txt",header = TRUE) %>% as.matrix
all91envir<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/all10pop_91environ_PCA_matrix.txt",header = TRUE) %>% as.matrix
Fst <- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/reorder_all10qa_fst_scalar_matrix.txt",header = TRUE) %>% as.matrix
fdpgi10kb<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/fd/fdpgi_all10qv_10k.txt",header = TRUE) %>% as.matrix
longitude<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/longitude_matrix.txt",header = TRUE) %>% as.matrix
latitude<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/latitude_matrix.txt",header = TRUE) %>% as.matrix
pre<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/all10pop_pre_PCA_matrix.txt",header = TRUE) %>% as.matrix
tem<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/all10pop_tem_PCA_matrix.txt",header = TRUE) %>% as.matrix
srad<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/all10pop_srad_PCA_matrix.txt",header = TRUE) %>% as.matrix
wind<- read.table("D:/麻栎数据/ABBA-BABA/方法二/equel2indsample/refilter_fd/matrix/PCA/all10pop_wind_PCA_matrix.txt",header = TRUE) %>% as.matrix


mantel(all91envir, fdpgi10kb)
mantel(geodistancematrix, fdpgi10kb)
mantel(longitude, fdpgi10kb)
mantel(latitude, fdpgi10kb)
mantel(pre, fdpgi10kb)
mantel(tem, fdpgi10kb)
mantel(srad, fdpgi10kb)
mantel(wind, fdpgi10kb)
mantel(all91envir, geodistancematrix)
# partial Mantel test
mantel.partial( geodistancematrix,fdpgi10kb, Fst)
mantel.partial( geodistancematrix,fdpgi10kb, all91envir)
mantel.partial(environ8matrix, fdpgi10kb, geodistancematrix)
mantel.partial(environ8matrix, fdpgi10kb, Fst)
mantel.partial(all91envir, fdpgi10kb, geodistancematrix)
mantel.partial(all91envir, fdpgi10kb, Fst)
mantel.partial(longitude, fdpgi10kb, geodistancematrix)
mantel.partial(longitude, fdpgi10kb, Fst)
mantel.partial(longitude, fdpgi10kb, all91envir)
mantel.partial(latitude, fdpgi10kb, geodistancematrix)
mantel.partial(latitude, fdpgi10kb, Fst)
mantel.partial(pre, fdpgi10kb, geodistancematrix)
mantel.partial(pre, fdpgi10kb, Fst)
mantel.partial(tem, fdpgi10kb, geodistancematrix)
mantel.partial(tem, fdpgi10kb, Fst)
mantel.partial(srad, fdpgi10kb, geodistancematrix)
mantel.partial(srad, fdpgi10kb, Fst)
mantel.partial(wind, fdpgi10kb, geodistancematrix)
mantel.partial(wind, fdpgi10kb, Fst)


#fst
mantel(geodistancematrix, Fst)
mantel(all91envir, Fst)
mantel(longitude, Fst)
mantel(latitude, Fst)
mantel(pre, Fst)
mantel(tem, Fst)
mantel(srad, Fst)
mantel(wind, Fst)
# partial Mantel test
mantel.partial( geodistancematrix,Fst, all91envir)
#mantel.partial(all91envir, Fst, geodistancematrix)
mantel.partial(all91envir, Fst, geodistancematrix)
#mantel.partial(longitude, Fst, geodistancematrix)
mantel.partial(longitude, Fst, geodistancematrix)
#mantel.partial(longitude, Fst, all91envir)
#mantel.partial(latitude, Fst, geodistancematrix)
mantel.partial(latitude, Fst, geodistancematrix)
#mantel.partial(pre, Fst, geodistancematrix)
mantel.partial(pre, Fst, geodistancematrix)
#mantel.partial(tem, Fst, geodistancematrix)
mantel.partial(tem, Fst, geodistancematrix)
#mantel.partial(srad, Fst, geodistancematrix)
mantel.partial(srad, Fst, geodistancematrix)
#mantel.partial(wind, Fst, geodistancematrix)
mantel.partial(wind, Fst, geodistancematrix)
