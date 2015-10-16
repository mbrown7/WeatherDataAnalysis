
##############################################################################
####FUNCTIONS
##############################################################################

#A function to subset the data into a single month-year combination
subsetData = function(month, year, data){
        monthU = month + 1
        if(month < 10){
                month = paste(0, month, sep = "")
        }
        if(monthU < 10){
                monthU = paste(0, monthU, sep = "")
        }
        lowerBound = paste(month, 0, 0, sep = "")
        upperBound = paste(monthU, 0, 0, sep = "")

        subset.df <- subset(data, data$YEAR == year)
        subset.df <- subset(subset.df, subset.df$MONTH_DAY > lowerBound)
        subset.df <- subset(subset.df, subset.df$MONTH_DAY < upperBound)
        subset.df
}

#A function to count the monthly occurrences of precipitation in different amount categories for a given year
monthlyPrecipitationCounts = function(year, data){
        dm <- matrix(0, 12, 6)
        row.names(dm) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        colnames(dm) <- c("0-0.25", "0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+")
        for(i in 1:12){
                subset.df <- subsetData(i, year, data)
                subset.hist.counts <- hist(subset.df$Precip.in, breaks = c(0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0), plot = FALSE)$counts
                dm[i, ] <- subset.hist.counts
        }
        dm
}

#A version of the above for temperature
monthlyTemp = function(year, data){
        dm <- matrix(0, 12, 15)
        row.names(dm) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
        colnames(dm) <- c("-300 to -250", "-250 to -200", "-200 to -150", "-150 to -100", "-100 to -50", "-50 to 0", "0 to 50",
	  	"50 to 100", "100 to 150", "150 to 200", "200 to 250", "250 to 300", "300 to 350", "350 to 400", "400 to 450")
        for(i in 1:12){
                subset.df <- subsetData(i, year, data)
                subset.hist.counts <- hist(subset.df$MaxTemp.0.1C, breaks = seq(-300, 450, 50), plot = FALSE)$counts
                dm[i, ] <- subset.hist.counts
        }
        dm
}

#A function to perform the 30-year-average, a climatological approach to missing data
fillMissingData = function(table, column){
	missingIndices <- which(is.na(table[, column]))

	for(i in 1:length(missingIndices)){
		year <- table[missingIndices[i], 2]
		md <- table[missingIndices[i], 3]
		#For a leap day, use Feb 28th
		if(md == "0229"){
			md = "0228"
		}
		#the earliest missing values we can fill are from 1901
		#if year is before the time when we have 30 year data
		if(year < 1901){
			cat("Error in Climatology: No 30-year Average\n")
			break
		}else if(year < 1911){
			#if year is between 1901 and 1910 inclusive, use 1871-1900
			lb <- 1871
			ub <- 1900
		}else if(year < 1921){
			#if year is between 1911 and 1920 inclusive, use 1881-1910
			lb <- 1881
			ub <- 1910
		}else if(year < 1931){
			lb <- 1891
			ub <- 1920
		}else if(year < 1941){
			lb <- 1901
			ub <- 1930
		}else if(year < 1951){
			lb <- 1911
			ub <- 1940
		}else if(year < 1961){
			lb <- 1921
			ub <- 1950
		}else if(year < 1971){
			lb <- 1931
			ub <- 1960
		}else if(year < 1981){
			lb <- 1941
			ub <- 1970
		}else if(year < 1991){
			lb <- 1951
			ub <- 1980
		}else if(year < 2001){
			lb <- 1961
			ub <- 1990
		}else if(year < 2011){
			lb <- 1971
			ub <- 2000
		}else{
			lb <- 2001
			ub <- 2010
		}

		table.s <- table[table$YEAR >= lb, ]
		table.s <- table.s[table.s$YEAR <= ub, ]

		s <- table.s[table.s$MONTH_DAY == md, ]

		table[missingIndices[i], column] <- mean(s[, column])
	}
	return(table)
}

#A function to portion the data into blocks of specified size (e.g. 10 years, 30 years, etc.)
makeIntoBlocks <- function(span, array, categories){
	x <- seq(minYear, maxYear, span)
	blocks <- array(data = NA, dim = c(12, categories, length(x)))
	for(i in 1:(length(x)-1)){
		a <- 1 + span * (i - 1)
		b <- span * i
		blocks[ , , i] <- apply(array[ , , a:b], MARGIN = c(1, 2), sum)
	}
	#the last case:
	a <- 1 + span * (length(x) - 1)
	b <- length(array[1, , ]) / categories
	blocks[ , , length(x)] <- apply(array[ , , a:b], MARGIN = c(1, 2), sum)
	blocks
}

#The same as above for Seasons
makeIntoBlocksSeasons <- function(span, array, categories){
	x <- seq(minYear, maxYear, span)
	blocks <- array(data = NA, dim = c(4, categories, length(x)))
	for(i in 1:(length(x)-1)){
		a <- 1 + span * (i - 1)
		b <- span * i
		blocks[ , , i] <- apply(array[ , , a:b], MARGIN = c(1, 2), sum)
	}
	#the last case:
	a <- 1 + span * (length(x) - 1)
	b <- length(array[1, , ]) / categories
	blocks[ , , length(x)] <- apply(array[ , , a:b], MARGIN = c(1, 2), sum)
	blocks
}

#Perform a proportion test for a month-year combination in the data set
	#x is the block category start years, such as 1869, 1899, etc.
	#set is the data frame we will use
	#month is the month in question (numeric: 1, 2, 3, ...)
	#y1 is the block designation (numeric: 1, 2, 3, ...), represents the number block (first, second, third, ...)
#The return value is a vector of p-values for proportion tests
#Each proportion test compares the data for a given month and a given block to the data for that month
#in all other blocks. This function needs to be called sequentially because it only compares to blocks
#that occur after the given block. (For instance, if y1 = 3 and there are 4 blocks, then the return
#will only compare blocks 3 and 4.)
propTests <- function(x, set, month, y1){
	p <- vector( )
	c <- 1
	a <- set[month, , y1]
	for(j in (y1 + 1):length(x)){
		b <- set[month, , j]
		for(k in 1:length(a)){
			p[c] <- prop.test(c(a[k], b[k]), c(sum(a), sum(b)))$p.value
			c <- c + 1
		}
	}
	p
}

#Prints out the significant blocks and categories from the propotion test results
#block is the number block that we are testing on
printSig <- function(w, block){ #made specifically for the size 10 and specifically for precipitation categories
	b <- ((w - 1) %/% 6) + block + 1 #the blocks
	m <- w %% 6
	for(i in 1:length(w)){
		y <- 1869 + 10 * (b[i] - 1)
		if(b[i] == 15){
			q <- maxYear
		}else{
			q <- y + 9
		}
		#z <- 0 + 50 * (m[i] - 1)
		#z is a bit different in the precipitation case since categories are not of equal size
		if(m[i] == 1){
			z1 <- 0
			z2 <- 0.25
		}else if(m[i] == 2){
			z1 <- 0.26
			z2 <- 0.5
		}else if(m[i] == 3){
			z1 <- 0.51
			z2 <- 1.0
		}else if(m[i] == 4){
			z1 <- 1.01
			z2 <- 2.0
		}else if(m[i] == 5){
			z1 <- 2.01
			z2 <- 3.0
		}else{ #in our heads m[i] is 6 but actually it would be 0, either way it works for this else
			z1 <- 3.0
			z2 <- "infinity"
		}
		cat("Block:", y, "to", q, "Section:", z1, "to", z2, "\n")
	}
}

#The same as above
printSig30 <- function(w, block){ #made specifically for the size 30
	b <- ((w - 1) %/% 6) + block + 1 #the blocks
	m <- w %% 6
	for(i in 1:length(w)){
		y <- 1869 + 30 * (b[i] - 1)
		if(b[i] == 5){
			q <- maxYear
		}else{
			q <- y + 29
		}
		if(m[i] == 1){
			z1 <- 0
			z2 <- 0.25
		}else if(m[i] == 2){
			z1 <- 0.26
			z2 <- 0.5
		}else if(m[i] == 3){
			z1 <- 0.51
			z2 <- 1.0
		}else if(m[i] == 4){
			z1 <- 1.01
			z2 <- 2.0
		}else if(m[i] == 5){
			z1 <- 2.01
			z2 <- 3.0
		}else{ #in our heads m[i] is 6 but actually it would be 0, either way it works for this else
			z1 <- 3.0
			z2 <- "infinity"
		}
		cat("Block:", y, "to", q, "Section:", z1, "to", z2, "\n")
	}
}

##############################################################################
####PLOTTING FUNCTIONS
##############################################################################

plotMonthlyEvents <- function(month, data, label){
	#Graph number of events in each category for each year in January
	d <- dimnames(data)[[3]]
	c <- colSums(data[month, , ])
	a <- matrix(0, 6, length(minYear:maxYear))
	for(i in 1:length(minYear:maxYear)){
		a[ , i] <- data[month, , i] / c[i]
	}
	monthName <- months[month]

	plot(d, a[1, ], xaxt = 'n', type = 'l', ylim = c(0, 1), main = paste(monthName, "Precipitation Events"),
		xlab = 'Year', ylab = label) #plot month in the 0-0.25 category witout axis
	points(d, a[2, ], xaxt = 'n', type = 'l', col = 'red') #month in the next precip category
	points(d, a[3, ], xaxt = 'n', type = 'l', col = 'green')
	points(d, a[4, ], xaxt = 'n', type = 'l', col = 'blue')
	points(d, a[5, ], xaxt = 'n', type = 'l', col = 'purple')
	points(d, a[6, ], xaxt = 'n', type = 'l', col = 'orange')
	axis(side = 1, at = c(seq(minYear, maxYear, 30), maxYear), labels = T) #add the X axis with a good number of tick marks
	legend("right", c("0-0.25", "0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+"),
		col = c('black', 'red', 'green', 'blue', 'purple', 'orange'), lwd = 2)
}

plotMonthlyEventsRainOnly <- function(month, data, label){
	#Graph number of events in each category for each year in January
	d <- dimnames(data)[[3]]
	c <- colSums(data[month, , ])
	a <- matrix(0, 6, length(minYear:maxYear))
	for(i in 1:length(minYear:maxYear)){
		a[ , i] <- data[month, , i] / c[i]
	}
	monthName <- months[month]
	plot(d, a[2, ], xaxt = 'n', type = 'l', ylim = c(0, max(a[2:6, ]) + 0.05), main = paste(monthName, "Precipitation Events"),
		xlab = 'Year', ylab = label) #plot January in the 0-0.25 category witout axis
	points(d, a[3, ], xaxt = 'n', type = 'l', col = 'green')
	points(d, a[4, ], xaxt = 'n', type = 'l', col = 'blue')
	points(d, a[5, ], xaxt = 'n', type = 'l', col = 'purple')
	points(d, a[6, ], xaxt = 'n', type = 'l', col = 'orange')
	axis(side = 1, at = c(seq(minYear, maxYear, 30), maxYear), labels = T) #add the X axis with a good number of tick marks
	legend("topright", c("0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+"),
		col = c('black', 'green', 'blue', 'purple', 'orange'), lwd = 2)
}

plotMonthlyLines <- function(month, data, label){
	#Graph number of events in each category for each year in January
	d <- as.numeric(dimnames(data)[[3]]) #all years
	c <- colSums(data[month, , ])
	a <- matrix(0, 6, length(minYear:maxYear))
	for(i in 1:length(minYear:maxYear)){
		a[ , i] <- data[month, , i] / c[i]
	}
	monthName <- months[month]
	plot(1869, 0, xaxt = 'n', type = 'l', xlim = c(minYear, maxYear), ylim = c(0, 1), main = paste(monthName, "Precipitation Events"),
		xlab = 'Year', ylab = label) #plot January in the 0-0.25 category witout axis
	axis(side = 1, at = c(seq(minYear, maxYear, 30), maxYear), labels = T) #add the X axis with a good number of tick marks

	for(i in 1:6){
		l <- lm(a[i, ] ~ d)
		abline(l$coeff[1], l$coeff[2], col = i)
	}
	legend("topright", c("0-0.25", "0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+"), col = c(1:6), lwd = 2)
}

plotMonthlyLinesRainOnly <- function(month, data, label){
	#Graph number of events in each category for each year in January
	d <- as.numeric(dimnames(data)[[3]]) #all years
	c <- colSums(data[month, , ])
	a <- matrix(0, 6, length(minYear:maxYear))
	for(i in 1:length(minYear:maxYear)){
		a[ , i] <- data[month, , i] / c[i]
	}
	monthName <- months[month]
	plot(1869, 0, xaxt = 'n', type = 'l', xlim = c(minYear, maxYear), ylim = c(0, 0.1), main = paste(monthName, "Precipitation Events"),
		xlab = 'Year', ylab = label) #plot January in the 0-0.25 category witout axis
	axis(side = 1, at = c(seq(minYear, maxYear, 30), maxYear), labels = T) #add the X axis with a good number of tick marks

	for(i in 2:6){
		l <- lm(a[i, ] ~ d)
		abline(l$coeff[1], l$coeff[2], col = i)
	}
	legend("right", c("0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+"), col = c(2:6), lwd = 2)
}

##############################################################################
####PRELIMINARIES and FORMATTING
##############################################################################

#Read the files
library(gdata)

s <- read.xls("C:/Users/Morgan/Desktop/School/Old/Summer15/IndepStudy/MadWeatherData.xlsx", 1) #fill the quotation marks with correct file path

#Set the minimum and maximum year here:
minYear = 1869
maxYear = 2014 #there is data for 2015 but the last complete year is 2014

s.df <- data.frame(s)

#Look for missing data
min(s.df$Precip.0.1mm) #This should be 0 if there is no missing data, but it is not 0
missingIndices <- which(s.df$Precip.0.1mm == -9999)
s.df[missingIndices, 3] = NA #Change the missing data locations to be NA

#Repeat for all other columns
s.df[which(s.df$SnowDepth.mm == -9999), 4] = NA
s.df[which(s.df$Snowfall.mm == -9999), 5] = NA
s.df[which(s.df$MaxTemp.0.1C == -9999), 6] = NA
s.df[which(s.df$MinTemp.0.1C == -9999), 7] = NA
s.df[which(s.df$WaterEquivalentSnowonGround.0.1mm == -9999), 8] = NA
s.df[which(s.df$PeakWindSpeedGust.0.1mps == -9999), 9] = NA

#Convert precipitation to inches
s.df$Precip.0.1mm <- s.df$Precip.0.1mm * 0.1
s.df$Precip.0.1mm <- s.df$Precip.0.1mm / 25.4
colnames(s.df)[3] <- "Precip.in"

#Separate the dates into year and month/day
s.df.sep <- data.frame(matrix(ncol = 10, nrow = length(s.df$DATE)))
colnames(s.df.sep) <- c("STATION_NAME", "YEAR", "MONTH_DAY", "Precip.in", "SnowDepth.mm", "Snowfall.mm",
	"MaxTemp.0.1C", "MinTemp.0.1C", "WaterEquivalentSnowonGround.0.1mm", "PeakWindSpeedGust.0.1mps")
s.df.sep[1] <- s.df[1]
s.df.sep[4] <- s.df[3]
s.df.sep[5] <- s.df[4]
s.df.sep[6] <- s.df[5]
s.df.sep[7] <- s.df[6]
s.df.sep[8] <- s.df[7]
s.df.sep[9] <- s.df[8]
s.df.sep[10] <- s.df[9]

x <- as.character(s.df$DATE)
y <- strsplit(x, "(?<=.{4})", perl = TRUE)
z <- lapply(y, "[[", 1)
s.df.sep$YEAR <- as.character(z)
z2 <- lapply(y, "[[", 2)
s.df.sep$MONTH_DAY <- as.character(z2)

#Now to see the entries which are NA you have to do
which(is.na(s.df.sep$Precip.in))

#Fill the missing data
s.df.sep <- fillMissingData(s.df.sep, 4) #the 4 argument specifies which column we want to fill

#Check that it worked
which(is.na(s.df.sep$Precip.in))

##############################################################################
####ORGANIZE the DATA
##############################################################################

#Create an array for data storage
allYears.precip.array <- array(data = NA, dim = c(12, 6, length(minYear:maxYear)), dimnames = list(c("Jan", "Feb", "Mar",
        "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), c("0-0.25", "0.26-0.5", "0.51-1.0", "1.01-2.0",
        "2.01-3.0", "3.0+"), c(minYear:maxYear)))

#Fill the array
for(i in 1:length(minYear:maxYear)){
        allYears.precip.array[ , , i] <- monthlyPrecipitationCounts(minYear - 1 + i, s.df.sep)
} #Because we are passing the whole object s.df.sep a bunch of times, this function takes very long to run

#To find out number of events in each category and each month across all years:
allYears.precip.array.aggregated <- apply(allYears.precip.array, MARGIN = c(1, 2), sum)

#Create an array for data storage based on the seasons
allYears.precip.array.seasonal <- array(data = NA, dim = c(4, 6, length(minYear:maxYear)), dimnames = list(c("Winter", "Spring",
        "Summer", "Fall"), c("0-0.25", "0.26-0.5", "0.51-1.0", "1.01-2.0", "2.01-3.0", "3.0+"), c(minYear:maxYear)))

#Fill the array
for(i in 1:length(minYear:maxYear)){
        allYears.precip.array.seasonal[2, , i] <- colSums(allYears.precip.array[3:5, , i]) #Spring
        allYears.precip.array.seasonal[3, , i] <- colSums(allYears.precip.array[6:8, , i]) #Summer
        allYears.precip.array.seasonal[4, , i] <- colSums(allYears.precip.array[9:11, , i]) #Fall
        allYears.precip.array.seasonal[1, , i] <- colSums(allYears.precip.array[1:2, , i]) + allYears.precip.array[12, , 1] #Winter
}


##############################################################################
####FIND SIGNIFICANT CHANGES
##############################################################################

#Note: ignore warnings from propTests( ) function

##########BLOCK SIZE 30

size <- 30
x <- seq(minYear, maxYear, size)
allYears.precip.blocks30 <- makeIntoBlocks(size, allYears.precip.array, 6) #There are 6 precipitation categories
pt <- propTests(x, allYears.precip.blocks30, 1, 1)
w <- which(pt < 0.05)
printSig30(w, 1)
#January

pt <- propTests(x, allYears.precip.blocks30, 2, 1)
w <- which(pt < 0.05)
printSig30(w, 1)
#February

pt <- propTests(x, allYears.precip.blocks30, 3, 1)
w <- which(pt < 0.05)
printSig30(w, 1)
#March

#...
#Use these to compare month-block combinations


##########SEASONAL, BLOCK SIZE 30

size <- 30
x <- seq(minYear, maxYear, size)
allYears.precip.seasonal.blocks30 <- makeIntoBlocksSeasons(size, allYears.precip.array.seasonal, 6) #There are 6 precipitation categories
pt <- propTests(x, allYears.precip.seasonal.blocks30, 1, 1)
w <- which(pt < 0.05)
printSig30(w, 1)
#Winter

pt <- propTests(x, allYears.precip.seasonal.blocks30, 2, 1)
w <- which(pt < 0.05)
printSig30(w, 1)
#Spring

#...

##########SEASONAL, BLOCK SIZE 10

size <- 10
x <- seq(minYear, maxYear, size)
allYears.precip.seasonal.blocks10 <- makeIntoBlocksSeasons(size, allYears.precip.array.seasonal, 6) #There are 6 precipitation categories
pt <- propTests(x, allYears.precip.seasonal.blocks10, 1, 1)
w <- which(pt < 0.05)
printSig(w, 1)

#...

##########SEASONAL, BLOCK SIZE 30

size <- 30
x <- seq(minYear, maxYear, size)
set <- makeIntoBlocksSeasons(size, allYears.precip.array.seasonal, 6)
ptt <- propTests(x, set, 1, 1)
w <- which(ptt < 0.05)
printSig30(w, 1)

#...

##############################################################################
####USING THE DATA FRAMES and NOTES
##############################################################################

#can use colSums to get the data for the whole year, like
colSums(allYears.precip.array[, , "1910"])
#gives the frequency of rainfall in different categories throughout all of 1910
#it never rained more than 3 inches in all of 1910

#note that rowSums when specifying year doesn't make much sense because
#it just gives the number of days in each month

#can also look at the occurences in a certain rain category for all years and months
allYears.precip.array[, "3.0+", ]

#can see the occurences just by year
colSums(allYears.precip.array[, "3.0+", ])
#i.e. in 1878 it rained more than 3.0 inches twice

#can see the occurences by month
rowSums(allYears.precip.array[, "3.0+", ])
#i.e. in the 42 year span it rained more than 3.0 inches in the warmer months (May through Sep) but not at all in any other
#month

#can look at data for a single month across all years
allYears.precip.array["Jun", , ]

#etc.... there are a lot of ways to look at the data and their relationships

#to obtain the relative frequency:
allYears.precip.array.aggregated.relative <- allYears.precip.array.aggregated / rowSums(allYears.precip.array.aggregated)

which.max(rowSums(allYears.precip.array.aggregated[, 2:6]))
#A rain event is defined as more than 0.25 inches
#the above command shows that the most rain events occur in June

#number of events in Spring
sum(allYears.precip.array.aggregated[3:5, 2:6])
#number of events in Summer
sum(allYears.precip.array.aggregated[6:8, 2:6])
#number of events in Fall
sum(allYears.precip.array.aggregated[9:11, 2:6])
#number of events in Winter
sum(allYears.precip.array.aggregated[c(1,2,12), 2:6])
#the most rain occurs in Summer, then Spring, then Fall, and then Winter

#The probability of rain events per month using 42 years of data
allYears.precip.array.aggregated.relative

#The above is interpreted as "the probability of it raining X inches given that it is Y month"

#NOTE the column names are a bit misleading because something that's 0.2567 would fall into the "0.26-0.5" category
#since it's more than the first break point of 0.25

##############################################################################
####VISUALIZING and WRITING the DATA
##############################################################################

hist(s.df.sep$Precip.in, breaks = c(0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)) #the graph is ugly
hist(s.df.sep$Precip.in, breaks = c(0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0))$counts

library(xlsx)
write.xlsx(allYears.precip.array.aggregated, "...") #write the data to file name "..." if desired

#Create a months object for help in naming files
months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
	"November", "December")

#Plot and save portion of events occurring in each rainfall category for each month
for(i in 1:12){
	png(filename = paste("...", months[i], ".png", sep = '')) #path name "..." is the folder where the file should be saved
	plotMonthlyEvents(i, allYears.precip.array, 'Portion of Events')
	dev.off( )
}

#Plot and save portion of events occurring in each rainfall category for each month
#Do not plot the 0-0.25 category
for(i in 1:12){
	png(filename = paste("...", months[i], ".png", sep = ''))
	plotMonthlyEventsRainOnly(i, allYears.precip.array, 'Portion of Events')
	dev.off( )
}

#Same as above but using only slope to indicate change in proportions over time
for(i in 1:12){
	png(filename = paste("...", months[i], ".png", sep = ''))
	plotMonthlyLines(i, allYears.precip.array, 'Portion of Events')
	dev.off( )
}
for(i in 1:12){
	png(filename = paste("...", months[i], ".png", sep = ''))
	plotMonthlyLinesRainOnly(i, allYears.precip.array, 'Portion of Events')
	dev.off( )
}