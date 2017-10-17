

#### OMXC20

data_stock <- read.table('data/omxc20.csv', header = T, sep = ';', dec = ',')
dim(data_stock)
head(data_stock)
str(data_stock)

plot(1:dim(data_stock)[1], data_stock$Closing.price)
hist(data_stock$Closing.price)


#### S&P 500