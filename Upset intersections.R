library(dplyr)

#convert any number other than 0 to 1
matrixupset[matrixupset != 0] <- 1

#already having the binary numbers and with unique columns, We will create a new df to work with.
newdf <- as.data.frame(matrixupset[, 5]) #create a new df with only the desired column
row.names(newdf) <- row.names(matrixupset) #that the rows are equal
colnames(newdf) <- colnames(matrixupset)[5] #that the columns are equal
newdf <- newdf %>% filter(newdf[,1] == 1) #filter only columns with 1

#We can extract unique proteins from columns with the code
unique_proteins <- matrixupset[rowSums(matrixupset[, -1]) == 0 & matrixupset$Control == 1, ]

#proteins in 2 variables, with the dplyr package
test1 <- matrixupset %>%
  filter(`Condition1` == 1 & `condition2` == 1) %>%
  filter(rowSums(select(., -c(2, 3, which(names(.) %in% c("Condition1", "condition2"))))) == 0)

#with 3 and so on
ControlSwitch1h4h <- matrixupset %>%
  filter(`Control` == 1 & `Condition1` == 1 & `condition2` == 1) %>%
  filter(rowSums(select(., -c(1, 2, 3, which(names(.) %in% c("Control", "Condition1", "condition2"))))) == 0)


####Take the intersections already having binary data of 1 for existence and 0 for not, THIS IS FOR 3 REPLICAS!!----
#Intersection in which the gene exists in all columns, only 2 of the replicas are needed.

tp10allin <- tp10[apply(tp10[, 1:3], 1, function(x) sum(x == 1) >= 2) & 
                    apply(tp10[, 4:6], 1, function(x) sum(x == 1) >= 2), ]

#intersection in which the gene exists in only 1 variable of the experiment
#in this case in the first 3 columns, where at least 2 are 1, and in the other 3 where at least 2 are 0.

tp9Ocontrol <- tp9[apply(tp9[, 1:3], 1, function(x) sum(x == 1) >= 2) & 
                     apply(tp9[, 4:6], 1, function(x) sum(x == 0) >= 2), ]

#at least 2 zeros in all 3 columns, and at least two 1s in the other 3 columns.

tp9starv <- tp9[apply(tp9[, 1:3], 1, function(x) sum(x == 0) >= 2) & 
                  apply(tp9[, 4:6], 1, function(x) sum(x == 1) >= 2), ]

#To extract row names in .txt

write(row.names(tp9Ocontrol), "tp9Ocontrol.txt")
write(row.names(tp9starv), "tp9starv.txt")

#in the txt you have to add "" to the genes so that Netlogo can read them
#is a bash to apply in Linux

for file in *.txt; do
temp_file=$(mktemp)
awk '{print "\"" $0 "\""}' "$file" | tr -d '\r' > "$temp_file"
mv "$temp_file" "$file"
done
