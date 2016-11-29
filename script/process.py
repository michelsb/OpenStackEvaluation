#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv     # imports the csv module
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
import numpy as np
#from rpy2.robjects.lib import ggplot2
import sys  

reload(sys)  
sys.setdefaultencoding('utf8')

types = ["ovs","bridge"]

rates = ["10000","30000","50000","70000","90000"]

r = robjects.r
r.library("nortest")
r.library("MASS") 

r('''
        wilcox.onesample.test <- function(v, verbose=FALSE) {
           wilcox.test(v,mu=median(v),conf.int=TRUE, conf.level = 0.99)
        }
        wilcox.twosamples.test <- function(v, r, verbose=FALSE) {
           wilcox.test(v,r, alternative = "greater")
        }
        find_outliers <- function(x, na.rm = TRUE, ...) {
  			qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  			H <- 1.5 * IQR(x, na.rm = na.rm)
  			y <- x
  			y[x < (qnt[1] - H)] <- NA
  			y[x > (qnt[2] + H)] <- NA
  			y
		}
		remove_na_values <- function(x) {  		
  			y <- x  		
  			y <- y[!is.na(y)]
  			y
		}
		remove_element <- function(x,index) {  		
  			y <- x[-index]  		
  			y
		}
        ''')

#r('''
#	remove_outliers <- function(x, na.rm = TRUE, ...) {  		
#  		y <- x  		
#  		y[x > 5000] <- NA
#  		y
#	}
#	''')

# Remove outliers
#find_outliers = robjects.r('find_outliers') 
#remove_na_values = robjects.r('remove_na_values')
#remove_element = robjects.r('remove_element')

# Normality Test
ks = robjects.r('ks.test') # KS
lillie = robjects.r('lillie.test') # Lilliefors
cvm = robjects.r('cvm.test') # Cramér-von Mises
shapiro = robjects.r('shapiro.test') # Shapiro-Wilk
sf = robjects.r('sf.test') # Shapiro-Francia
ad = robjects.r('ad.test') # Anderson-Darling

# Non-parametric Tes
wilcoxon_test_two_samples = robjects.r['wilcox.twosamples.test']#robjects.r('wilcox.test')
wilcoxon_test_one_sample = robjects.r['wilcox.onesample.test']

vectors_mean = dict()
vector_samples_per_rate = dict()

for typ in types:

	vectors_mean[typ] = []
	means = vectors_mean[typ]
	#vectors_median[typ] = []
	#medians = vectors_median[typ]	

	if typ == "ovs":
		title = "OpenVSwitch "
	else:
		title = "Linux Bridge "

	for rate in rates:

		if str(rate) not in vector_samples_per_rate.keys():
			vector_samples_per_rate[str(rate)] = dict()

		f = open('../results/stats-'+typ+'-2-512-'+rate+'.log', 'rb') # opens the csv file		

		vector = []

		try:
		    reader = csv.reader(f,delimiter=' ') #creates the reader object
		    for row in reader:   # iterates the rows of the file in orders		    	
		        vector.append(float(row[2])*1000000)    # prints each row
		finally:
		    f.close()      # closing

		v = robjects.FloatVector(vector)
		
		#v = find_outliers(v)

		if str(rate) not in vector_samples_per_rate.keys():
			vector_samples_per_rate[str(rate)] = dict()

		vector_samples_per_rate[str(rate)][typ]=v 

		#v = remove_na_values(v)
		
		title_graph = title + "(" + str(rate) + " pps)"

		r.pdf("../figures/hist-"+typ+"-2-512-"+rate+".pdf")
		r.hist(v, main = title_graph, col="blue", xlab = "Latência (microssegundos)", ylab = "Frequência Absoluta")

		r.pdf("../figures/box-"+typ+"-2-512-"+rate+".pdf")
		r.boxplot(v, main = title_graph,col="lightblue", horizontal=True, las=1, xlab="Latência (microssegundos)")

		# Grafico de probabilidade (QQ)
		r.pdf("../figures/qq-"+typ+"-2-512-"+rate+".pdf")
		r.qqnorm(v, main = title_graph, xlab = "Quantis teóricos N(0,1)", pch = 20,
		 ylab = "Latência (microssegundos)")
		r.qqline(v, lty = 2, col = "red")		

		mean = r.mean(v)[0]
		sd = r.sd(v)[0]		

		print "------- Statistics "+ typ +" "+ rate +" -------"		
		print r.summary(v)
		print "-------------------------- \n"

		means.append(r.mean(v)[0])
		#medians.append(r.median(v)[0])		

		norm_tests = open('../norm-tests/norm-tests-'+typ+'-2-512-'+rate+'.csv', 'wb') # opens the csv file		

		try:			
			writer = csv.writer(norm_tests, delimiter=' ', quoting=csv.QUOTE_MINIMAL)			
			writer.writerow(['Method', 'Statistic', 'P-Value','Alternative Hypothesis (KS Test)'])

			test = ks(v,"pnorm",mean,sd)			
			writer.writerow([test[3][0], test[0][0],test[1][0],test[2][0]])

			test = lillie(v)
			writer.writerow([test[2][0], test[0][0],test[1][0],''])
			test = cvm(v)
			writer.writerow([test[2][0], test[0][0],test[1][0],''])
			test = shapiro(v)
			writer.writerow([test[2][0], test[0][0],test[1][0],''])
			test = sf(v)
			writer.writerow([test[2][0], test[0][0],test[1][0],''])
			test = ad(v)
			writer.writerow([test[2][0], test[0][0],test[1][0],''])
		finally:
		    norm_tests.close()		
		


#plot_ovs, =  plt.plot( x, vectors_mean["ovs"], 'go', label='OpenVSwitch') # green bolinha
#plt.plot( x, vectors_mean["ovs"], 'k:', label='OpenVSwitch', color='orange') # linha pontilha orange

#plot_bridge, = plt.plot( x, vectors_mean["bridge"], 'r^', label='Linux Bridge') # red triangulo
#plt.plot( x, vectors_mean["bridge"], 'k--', label='Linux Bridge', color='blue')  # linha tracejada azul

#plt.legend(handles=[plot_bridge, plot_ovs])
#plt.grid(True)
#plt.xlabel("Taxa de Transmissão (pps)")
#plt.ylabel("Latência (microssegundos)")
#plt.savefig("../figures/mean-ovs-bridge-2-512.png")

vectors_median = {'ovs':[],'bridge':[]}
errors_median = {'ovs':[],'bridge':[]}

try:

	np_tests = open('../np-tests/np-tests-wilcoxon-2-512.csv', 'wb')	
	writer = csv.writer(np_tests, delimiter=' ', quoting=csv.QUOTE_MINIMAL)
	writer.writerow(['Rate', 'Statistic', 'P-Value','Alternative Hypothesis'])	

	for rate in rates:	
		
		#count = 1
		#for value in vector_samples_per_rate[str(rate)]['bridge']:
		#	if str(value) == "NA":
		#		vector_samples_per_rate[str(rate)]['bridge'] = remove_element(vector_samples_per_rate[str(rate)]['bridge'],count)
		#		vector_samples_per_rate[str(rate)]['ovs'] = remove_element(vector_samples_per_rate[str(rate)]['ovs'],count)
		#	else:
		#		count += 1

		#count = 1
		#for value in vector_samples_per_rate[str(rate)]['ovs']:
		#	if str(value) == "NA":
		#		vector_samples_per_rate[str(rate)]['bridge'] = remove_element(vector_samples_per_rate[str(rate)]['bridge'],count)
		#		vector_samples_per_rate[str(rate)]['ovs'] = remove_element(vector_samples_per_rate[str(rate)]['ovs'],count)
		#	else:
		#		count += 1

		test = wilcoxon_test_one_sample(vector_samples_per_rate[str(rate)]['bridge'])		
		error_max = test[7][1]		
		median = test[8][0]
		vectors_median['bridge'].append(median)		
		errors_median['bridge'].append(float(error_max)-float(median))
		test = wilcoxon_test_one_sample(vector_samples_per_rate[str(rate)]['ovs'])	
		error_max = test[7][1]		
		median = test[8][0]
		vectors_median['ovs'].append(median)	
		errors_median['ovs'].append(float(error_max)-float(median))
		np_test = wilcoxon_test_two_samples(vector_samples_per_rate[str(rate)]['bridge'],vector_samples_per_rate[str(rate)]['ovs'])#,paired = True)		
		writer.writerow([rate, np_test[0][0],np_test[2][0],np_test[4][0]])

finally:
	np_tests.close()

# Drawing median an mean graphics
x = np.array(rates)

plot_ovs, =  plt.plot( x, vectors_median["ovs"], 'go', label='OpenVSwitch') # green bolinha
plt.plot( x, vectors_median["ovs"], 'k:', label='OpenVSwitch', color='orange') # linha pontilha orange
plt.errorbar(x, vectors_median["ovs"], yerr=errors_median['ovs'], linestyle="None", marker="None", color="orange")

plot_bridge, = plt.plot( x, vectors_median["bridge"], 'r^', label='Linux Bridge') # red triangulo
plt.plot( x, vectors_median["bridge"], 'k--', label='Linux Bridge', color='blue')  # linha tracejada azul
plt.errorbar(x, vectors_median["bridge"], yerr=errors_median['bridge'], linestyle="None", marker="None", color="blue")

plt.title("Gráfico das Medianas")
plt.legend(handles=[plot_bridge, plot_ovs])
plt.grid(True)
plt.xlabel("Taxa de Transmissão (pps)")
plt.ylabel("Latência (microssegundos)")
plt.savefig("../figures/median-ovs-bridge-2-512.pdf",format='pdf')
