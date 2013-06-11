
### How to call this

## Either you let the run_calls.pl script do everything. It will call this, and pass the output to process_call.pl,
##            and you never need to think about this.

## OR, you run this yourself, and then call process_calls.pl

##   cat classifier.parallel.ploidy_aware.R | R --vanilla --args arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9


#### Prerequisite: 1. you need covg and read length info per colour. This is printed to out by Cortex when running 
####                  I assume you have saved this in a file. From this, you need to generate a little file containing a table of this data
###                   If the cortex output is in log.txt, then run :  perl make_read_len_and_total_seq_table.pl log.txt > table.txt
###                   This will be arg5 for this R script

### Prerequisite:  2. I don't want to have to parse a Cortex output file in R, so I have a script make_covg_file.pl 
###                   which makes an appropriate input file for the classifier.
###                   You do this:  

#     perl /cortex release dir/scripts/analyse_variants/make_covg_file.pl <cortex callfile> 
#                                                                         <number of colours in graph> 
#                                                                         <ref colour> (using -1 if no ref) > covg_for_classifier_file


### Prerequisite 3:   Since the classifier is a bit slow when running on millions of variants
###                   this is set up to allow you to run many copies of this in parallel, and then concatenate them into one file at the end
###                   You tell it on the commandline what is the first variant to start with, and how many to do, and give it an input file
###                   If you don't care, or if this is a small dataset, or if this is the first time you have done it, I suggest you just
###                   use num_vars_to_process (arg2) =totalvars(arg4) = however many variants are in your file.


### arg 1 - number of the first variant to use (if the first one is var_1, then enter 1).
### arg 2 - how many variants to process/classify
### arg 3 - input covg_for_classifier file. You should generate this in advance
### arg 4 - number of rows/lines in the covg_for_classifier file = number of variants overall
### arg 5 - number of colours in the graph 
###         For example  you use cortex_var_31_c7 to get your calls, then this argument should be 7, even if you only loaded data from one sample in
### arg 6 - was there a reference colour? 1 for yes and 0 for no (doesn't matter which colour it was)
### arg 7  - table of read lengths and covgs
### arg 8  - estimated genome size. (Don't panic if not exact)
### arg 9  - kmer size
### arg 10 - ploidy. 1 for haploid, 2 for diploid, no other value acceptable.
### arg 11 - output file name




args <- commandArgs(trailingOnly = TRUE) ##  one argument only -  start number 

first_var            <- as.integer(args[1]);
num_vars_to_process  <- as.integer(args[2])
covgfile             <- args[3]
totalvars            <- as.integer(args[4])
num_colours          <- as.integer(args[5])
ref_present          <- as.integer(args[6])
tablefile            <- args[7]
genome.size          <- as.integer(args[8]);
kmer                 <- as.integer(args[9]);
ploidy               <- as.integer(args[10]);
outfile              <- args[11];

## covg file will contain two char columns (var name, whether this looks like a ref bubble) which are 
## not used by the classifier, followed by 2 numeric columns (length branch1, length branch2)
## followed by 2*num_colours numeric columns (covg on the two branches in each colour)
## If there was a reference in the graph, then the first two of these are covgs in the ref, also not used by the classifier.

num_numeric_cols            <- 2*num_colours
first_column_of_sample_data <- 5;
num_samples                 <- num_colours

if (ref_present==1)
{
    ## ignore the two columns showing covg in ref - that's purely for humans wanting to debug/explore.
    ## first 6 columns are: var name, is it a ref bubble, branch1 lenght, branch2 length, ref covg br1, ref covg br2
	first_column_of_sample_data<- 7; 
	num_samples                <- num_colours-1
}


d<-read.table(covgfile, as.is=T, head=F,
 		nrows= totalvars, colClasses=c("character", "character", rep("numeric",2*num_colours + 2)), 
		comment.char = "", sep="\t"); 


#outfile<-paste(covgfile, "split_start_", first_var, sep="")
logfile              <- paste(outfile, ".log", sep="");
last_var             <- first_var + num_vars_to_process -1



## The following  are the numbers of base-pairs loaded into each colour - take these from your Cortex output
stats_table<-read.table(tablefile, header=FALSE)
cov <- as.numeric(unlist(stats_table[,3]))


#### read length
read.length<-as.numeric(unlist(stats_table[,2]));



#Function to simulate Poisson data with var = 2*mean using gamma for lambda

r.compound.pois<-function(n=1, mean=1, f=2) {
	beta<-1/(f-1);
	alpha<-mean*beta;
	lambda<-rgamma(n, alpha, beta);
	return(rpois(n, lambda));
}


#alpha.bal is the symmetric coefficient for the balance prior
#gamma.rtp is the geometric parameter for the rpt cpt in the genome.  Mn = 1/(gamma.rpt)+1



llk.cov<-function(data, a.cov=0, g.min=1, g.max=10, over.dispersed=TRUE, f=2, rel.cov=rep(1, nrow(data))) {

	if (f<=1) {
		cat("\n\nError: Cannot have inflation of variance < 1 - setting to 1\n\n", file=logfile, append=TRUE,  sep="\n");
		f<-1.001;
	}

	n<-nrow(data);

	#Estimate mean allele coverage if not known
	if (a.cov==0) {
		a.cov<-sum(data)/(2*n);
	}

	g.range<-g.min:g.max;
	
	alpha<-2*g.range*a.cov/(f-1);
	beta<-alpha;
	sigma<-sum(data);
	cts<-apply(data, 1, sum);

	#V1

	if (over.dispersed) {
		llk.cov.g<-sigma*log(g.range)+alpha*log(beta)-lgamma(alpha)+lgamma(alpha+sigma)-
			(alpha+sigma)*log(beta+2*g.range*a.cov*n);
		llk.cov.g<-llk.cov.g-max(llk.cov.g);

		return(cbind(g.range, llk.cov.g));
	}

	#V3 - Pure Poisson
	if (!over.dispersed) {
		llk.cov.g<-sigma*log(g.range)-2*g.range*a.cov*n;
		llk.cov.g<-llk.cov.g-max(llk.cov.g);

		return(cbind(g.range, llk.cov.g));
	}

}

#To calculate likelihood of allele balance under repeat model

llk.bal.rep<-function(data, alpha.bal=5) {
	cts<-apply(data, 2, sum);
	return(lbeta(alpha.bal+cts[1], alpha.bal+cts[2])-lbeta(alpha.bal, alpha.bal));
}

llk.bal.err<-function(data, err=0.05, scale=100) {
	cts<-sort(apply(data, 2, sum));
	vals<-c(scale*err, scale);
	return(lbeta(vals[1]+cts[1], vals[2]+cts[2])-lbeta(vals[1],vals[2]));
}

llk.bal.snp<-function(data, grid.step=1/(2*nrow(data)), err=0.05) {

	#Define prior for allele and genotype frequencies
	f.grid<-seq(grid.step, 1-grid.step, grid.step);
	p.f<-1/(f.grid*(1-f.grid));
	p.f<-p.f/sum(p.f);
	n<-nrow(data);
	p.gt<-cbind((1-f.grid)^2, 2*f.grid*(1-f.grid), f.grid^2);
	gs<-c(2,1,0);
	if (ploidy==1)
	{
		p.gt<-cbind(1-f.grid, f.grid);	
		gs<-c(2,0);
	}	

	es<-(gs*(1-err)+(2-gs)*err)/2;

	llk.a.1<-data[,1, drop=F] %*% log(es);
	llk.a.2<-data[,2, drop=F] %*% log(1-es);
	llk.a<-llk.a.1+llk.a.2;
	mx<-apply(llk.a, 1, max);
	lk.a<-t(exp(llk.a-mx));

	#Get LLks for each allele frequency
	llk.f<-rep(0, nrow(p.gt));
	for (i in 1:nrow(p.gt)) llk.f[i]<-sum(log(p.gt[i,,drop=F]%*%lk.a));
	llk.f<-llk.f+sum(mx);

	#Include prior and integrate	
	mx<-max(llk.f);
	llk<-mx+log(sum(p.f*exp(llk.f-mx)));

	return(llk);

}


llk.cov.rep<-function(input, p.rep=0.5) {
	ps<-dgeom(input[2:nrow(input),1]-1, p.rep);
	ps<-ps/sum(ps);
	mx<-max(input[2:nrow(input),2]);
	vals<-exp(input[2:nrow(input),2]-mx);
	return(log(sum(vals*ps))+mx);
}





#x11()

a.cov=cov/(genome.size);

k.cov <- rep(0, length(read.length))
is.zero <- (read.length == 0)
k.cov[!is.zero] <- a.cov[!is.zero]*(read.length[!is.zero]-kmer+1)/read.length[!is.zero];

per.base.read.arrival<-rep(0, length(read.length))
per.base.read.arrival[!is.zero]<-cov[!is.zero]/(read.length[!is.zero]*genome.size);

#First get coverage for each sample from bubbles - check looks concordant
cov.mn<-apply(d[,first_column_of_sample_data:ncol(d)], 2, mean, trim=0.05);

ii<-1:num_samples;

rel.cov<- rep(0, length(per.base.read.arrival))
rc.is.zero <- (rel.cov==0)
mn.cov<-mean(per.base.read.arrival);
rel.cov[!rc.is.zero]<-per.base.read.arrival[!rc.is.zero]/mn.cov;


cat("\n\n***Classifying bubbles***\n\n", file=logfile, append=TRUE,  sep="\n");

to.do<-nrow(d);
num_cols_in_op<-7   
op<-array(0, c(to.do, num_cols_in_op));
colnames(op)<-c("ID", "LLk.G1", "LLk.SNP", "LLk.ERR", "LLk.REP", "REF.BUBBLE", "Model");

lim<-min(last_var, to.do)


#for (i in 1:to.do) {
for (i in first_var:lim) {


        cat("\n", file=logfile, append=TRUE,  sep="");
	cat("Doing site ", i, file=logfile, append=TRUE,  sep="");
	d1<-t(array(as.integer(d[i,first_column_of_sample_data:ncol(d)]), c(2, num_samples)));

	k.cov.emp<-(d[i,3]+d[i,4])*mn.cov/2;

	llk1<-llk.cov(d1, a.cov=k.cov.emp, rel.cov=rel.cov, f=2);
	llk.cov.rep1<-llk.cov.rep(llk1);

	llk.bal.rep1<-llk.bal.rep(d1, alpha.bal=2);
	llk.bal.err1<-llk.bal.err(d1, err=0.01);
	llk.bal.snp1<-llk.bal.snp(d1, err=0.01);

	op[i,1]<-i
	op[i,2]<-llk1[1,2];

	op[i,3]<-op[i,2]+llk.bal.snp1;
	op[i,4]<-op[i,2]+llk.bal.err1;
	#op[i,5]<-llk.cov.rep1 + llk.bal.rep1; 
	op[i,5]<- max(llk1[,2]) + llk.bal.rep1; 
	op[i,6]<-as.integer(d[i,2]=="REF_BUBBLE");

	a1<- d1[,1]>0;
	a2<- d1[,2]>0;
	max.gt<-rep(0, nrow(d1));
	max.gt[a1 & a2]<-1;
	max.gt[!a1 & a2]<-2;
	gt.cts<-hist(max.gt, plot=F, breaks=c(-0.5,0.5,1.5,2.5))$counts


	op[i,7]<-which.max(op[i,3:5]);

	
}
cat("\n\nDone!\n\n", file=logfile, append=TRUE,  sep="\n");

thresh<-log(10);
#min.cov<-8;
llk.thresh<- -40;
#cv.cov.thresh<-10;

no.cov<-apply((d[1:to.do,2*ii+3, drop=FALSE]+d[1:to.do,2*ii+4, drop=FALSE])>0, 1, sum);
var.cov<-apply((d[1:to.do,2*ii+3, drop=FALSE]+d[1:to.do,2*ii+4, drop=FALSE]), 1, var);
cv.cov<-var.cov/apply((d[1:to.do,2*ii+3, drop=FALSE]+d[1:to.do,2*ii+4, drop=FALSE]), 1, mean);
#which.snp<-which((op[,7]==1) & ((op[,3]-apply(op[,4:5], 1, max))>thresh) & 
#	(no.cov>=min.cov) & op[,3]>llk.thresh & cv.cov<cv.cov.thresh & op[,7]>hwe.thresh);

which.snp<-which((op[,7]==1) & ((op[,3]-apply(op[,4:5], 1, max))>thresh) )
## for debugging comment out which.snp<-which((op[,7]==1) & ((op[,3]-apply(op[,4:5, DROP=FALSE], 1, max))>thresh) ) 


############################################################
### Optional - do some plots
############################################################

#x11()

if (0)
{
plot(cov.mn[2*ii-1], cov.mn[2*ii], pch=19, col="blue", xlim=c(0,max(cov.mn)), ylim=c(0,max(cov.mn)) );
abline(0,1,col="red", lty="dotted");  #Highly correlated
samp.cov<-cov.mn[2*ii-1]+cov.mn[2*ii];
}
if (0)
{
plot(k.cov, samp.cov, pch=19, col="blue", xlab="Predicted kmer cov", 
	ylab="Mean kmer cov at bubbles", xlim=c(0,max(k.cov)), ylim=c(0, max(samp.cov)));
abline(0,1,col="red", lty="dotted");
}

if (0)
{
#Hists of coverage
par(mfrow=c(2,5));
rng<-quantile(unlist(d[,3:ncol(d)]), 0.9);
for (i in 1:10) { #just plot 10 of the num_samples
	hist(d[,2*i+5]+d[,2*i+6], col="blue", xlab="kmer ct", xlim=c(0,10), breaks=c(0:10, 1000),
	main=paste("kCOV = ", signif(samp.cov[i], 2), sep=""));
}
#par(mfrow=c(1,1));
}
if (0)
{
#QQplots against overdispersed Poisson
par(mfrow=c(2,5));
for (i in 1:10) {
	qqplot(r.compound.pois(1000, samp.cov[i], f=2), d[,2*i+5]+d[,2*i+6], pch=19, 
		xlab="Compound Poisson", ylab="Data");
	abline(0,1,col="red", lty="dotted");
 }
}



#######################
### output useful stats
#######################
cat("\n", sep="", file=logfile, append=TRUE);
cat("\nClassify ", 100*length(which.snp)/num_vars_to_process, "% as variants", sep="", file=logfile, append=TRUE);
#cat("\nOf those passing threshold, ", 100*mean(op[which.snp,6]), 
#	"% are REF BUBBLE compared to ", 100*mean(op[,6]), "% overall\n\n", sep="", file=logfile, append=TRUE);
#which.fp<-which.snp[op[which.snp,6]==1];
#which.rep<-which(op[,7]==3);
#which.err<-which(op[,7]==2);
#cat("\n\nOf ", length(which.rep), " bubbles identified as repeat, ", 100*mean(op[which.rep,6]), 
#	"% are REF bubbles\n\n", sep="", file=logfile, append=TRUE);
#cat("\n\nOf ", length(which.err), " bubbles identified as errors, ", 100*mean(op[which.err,6]), 
#	"% are REF bubbles\n\n", sep="", file=logfile, append=TRUE);
#ss<-setdiff(which(op[,7]==1), which.snp);
#cat("\n\nOf ", length(ss), " bubbles identified as variants, but fail to meet threshold, ", 100*mean(op[ss,6]), 
#	"% are REF bubbles\n\n", sep="", file=logfile, append=TRUE);
#cat("\n\nClassification eliminates ", 100-100*sum(op[which.snp,6])/sum(op[,6]), "% of errors\n\n", sep="\n", file=logfile, append=TRUE);



##############################################
## output results to file
##############################################

## for each call, print variant name \t winning model \t confidence
#for (i in 1:to.do) 


for (i in first_var:lim) 
{
	cat(d[op[i,1],1], "\t", file=outfile, append=TRUE); 
	if (abs(op[i,7]-1)<0.5)
	{
		cat("variant\t", op[i,3]-max(op[i,4:5]),  file=outfile, append=TRUE,  sep=""); 
	}
	else if  (abs(op[i,7]-2)<0.5)
        {
               cat("error\t", op[i,4]-max(c(op[i,3],op[i,5])),	file=outfile, append=TRUE,  sep=""); 
	}
	else if (abs(op[i,7]-3)<0.5)
	{
		cat("repeat\t", op[i,5]-max(op[i,3:4]),  file=outfile, append=TRUE,  sep="");	 
	}
	else
	{
		cat("not classified", file=outfile,append=TRUE, sep=""); 
	}
	cat("\n", file=outfile, append=TRUE, sep="");
}



