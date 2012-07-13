## Search for ##<<<<<<<<<<<< MODIFY
##  tells you what you need to change


## modify this to hold the output of make_covg_file.pl   ##<<<<<<<<<<<< MODIFY
d<-read.table("cortex_callfile.covg_for_classifier.strip_header", as.is=T, head=F);


## The following  are the numbers of base-pairs loaded into each colour - take these from your Cortex output
## For example - suppose when you built your Cortex binary you got this output

#****************************************
#SUMMARY:
#Colour: MeanReadLen     TotalSeq
#0       0       2912760135
#1       50      18286122352
#2       50      16816361244
#3       50      18039181209
#4       50      15879192506
#5       50      17729089947
#6       50      15750659112
#7       50      26196361173
#8       52      20202087523
#9       50      18907785783
#10      50      16870486574
#****************************************

### then you fill in those numbers here - not I ignore the ref colour, 0
cov<-c(18286122352, 16816361244, 18039181209, 15879192506, 17729089947, 15750659112, 
	26196361173, 20202087523, 18907785783, 16870486574);##<<<<<<<<<<<< MODIFY

#genome length##<<<<<<<<<<<< MODIFY
genome.size<-3e9;


#### read length<<<<<<<<<<<< MODIFY
read.length<-50;

## kmer size
kmer<-31; ##<<<<<<<<<<<< MODIFY


#number of samples  ##<<<<<<<<<<<< MODIFY
num_samples<-10;

#Function to simulate Poisson data with var = 2*mean using gamma for lambda

r.compound.pois<-function(n=1, mean=1, f=2) {
	beta<-1/(f-1);
	alpha<-mean*beta;
	lambda<-rgamma(n, alpha, beta);
	return(rpois(n, lambda));
}

#Function to calculate integrated likelihood of REP model given data
#alpha.bal is the symmetric coefficient for the balance prior
#gamma.rtp is the geometric parameter for the rpt cpt in the genome.  Mn = 1/(gamma.rpt)+1
#NB should really redo with incomplete gamma function


llk.rep<-function(data, err=0.01, a.cov=0, alpha.bal=5, gamma.rpt=0.5, g.max=10, 
	over.dispersed=TRUE, f=2) {

	n<-nrow(data);

	#Estimate mean allele coverage if not known
	if (a.cov==0) {
		a.cov<-sum(data)/(2*n);
	}
	c1<-sum(data[,1]);
	c2<-sum(data[,2]);
	c<-c1+c2;

	#LLk from allele balance
	llk.bal<- -lbeta(alpha.bal, alpha.bal)+lbeta(alpha.bal+c1, alpha.bal+c2);
	llk.bal<-llk.bal +sum(lchoose(data[,1]+data[,2], data[,1]));

	#LLk from coverage
	p.g<-dgeom(0:(g.max-1), gamma.rpt);
	p.g<-p.g/sum(p.g);

	if (!over.dispersed) {
			llk.cov.g<-lapply(data[,1]+data[,2], dpois, a.cov*2*(1:g.max), log=T);
	}
	llk.cov.g<-array(unlist(llk.cov.g), c(n, g.max));
	llk.cov.g<-apply(llk.cov.g, 1, sum);
	mx<-max(llk.cov.g);
	sm<-sum(p.g*exp(llk.cov.g-mx));
	llk.cov<-log(sm)+mx;

	return(llk.bal+llk.cov);

}


llk.cov<-function(data, a.cov=0, g.min=1, g.max=10, over.dispersed=TRUE, f=2, rel.cov=rep(1, nrow(data))) {

	if (f<=1) {
		cat("\n\nError: Cannot have inflation of variance < 1 - setting to 1\n\n");
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
	p.gt<-cbind((1-f.grid)^2, 2*f.grid*(1-f.grid), f.grid^2);
	n<-nrow(data);

	gs<-c(2,1,0);
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






#Function to calculate integrated likelihood of SNP model given data
#a.cov can be estimated or inputted

llk.snp<-function(data, err=0.01, a.cov=0, grid.step=0.05) {

	#Define prior for allele and genotype frequencies
	f.grid<-seq(grid.step, 1-grid.step, grid.step);
	p.f<-1/(f.grid*(1-f.grid));
	p.f<-p.f/sum(p.f);
	p.gt<-cbind((1-f.grid)^2, 2*f.grid*(1-f.grid), f.grid^2);
	n<-nrow(data);

	#Estimate mean allele coverage if not known
	if (a.cov==0) {
		a.cov<-sum(data)/(2*n);
	}
	c1<-sum(data[,1]);
	c2<-sum(data[,2]);
	c<-c1+c2;

	#Genotypes and expected allele counts given genotypes allowing for error
	gs<-c(2,1,0);
	es<-gs*(1-err)+(2-gs)*err;

	#Get LLks for each allele in each sample and combine for a per sample per GT LLk
	llk.a.1<-lapply(data[,1], dpois, es*a.cov, log=T);
	llk.a.1<-array(unlist(llk.a.1), c(3, n));
	llk.a.2<-lapply(data[,2], dpois, (2-es)*a.cov, log=T);
	llk.a.2<-array(unlist(llk.a.2), c(3, n));
	llk.a<-t(llk.a.1+llk.a.2);
	mx<-apply(llk.a, 1, max);
	lk.a<-t(exp(llk.a-mx));

	#Get LLks for each allele frequency
	llk.f<-rep(0, nrow(p.gt));
	for (i in 1:nrow(p.gt)) llk.f[i]<-sum(log(p.gt[i,]%*%lk.a));
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

a.cov<-cov/(genome.size);
k.cov<-a.cov*(read.length-kmer+1)/read.length;

per.base.read.arrival<-cov/(read.length*genome.size);

#First get coverage for each sample from bubbles - check looks concordant
cov.mn<-apply(d[,7:ncol(d)], 2, mean, trim=0.05);

ii<-1:num_samples;
#plot(cov.mn[2*ii-1], cov.mn[2*ii], pch=19, col="blue");
#abline(0,1,col="red", lty="dotted");  #Highly correlated
samp.cov<-cov.mn[2*ii-1]+cov.mn[2*ii];
#plot(k.cov, samp.cov, pch=19, col="blue", xlab="Predicted kmer cov", 
#	ylab="Mean kmer cov at bubbles");
#abline(0,1,col="red", lty="dotted");

#Hists of coverage
#par(mfrow=c(2,5));
rng<-quantile(unlist(d[,3:ncol(d)]), 0.9);
#for (i in 1:10) {
#	hist(d[,2*i+5]+d[,2*i+6], col="blue", xlab="kmer ct", xlim=c(0,10), breaks=c(0:10, 1000),
#	main=paste("kCOV = ", signif(samp.cov[i], 2), sep=""));
#}
#par(mfrow=c(1,1));

#QQplots against overdispersed Poisson
#par(mfrow=c(2,5));
#for (i in 1:10) {
#	qqplot(r.compound.pois(1000, samp.cov[i], f=2), d[,2*i+5]+d[,2*i+6], pch=19, 
#		xlab="Compound Poisson", ylab="Data");
#	abline(0,1,col="red", lty="dotted");
#}
#par(mfrow=c(1,1));

rel.cov<-per.base.read.arrival/mean(per.base.read.arrival);
mn.cov<-mean(per.base.read.arrival);

cat("\n\n***Classifying bubbles***\n\n");

#Just look at SNPS
#d<-d[d[,3]==32 & d[,4]==32,];

to.do<-nrow(d);
op<-array(0, c(to.do, 8));
colnames(op)<-c("ID", "LLk.G1", "LLk.SNP", "LLk.ERR", "LLk.REP", "REF.BUBBLE", "p.HWE", "Model");

for (i in 1:to.do) {

	cat("\rDoing site ", i);
	d1<-t(array(as.integer(d[i,7:ncol(d)]), c(2, 10)));

	k.cov.emp<-(d[i,3]+d[i,4])*mn.cov/2;

	llk1<-llk.cov(d1, a.cov=k.cov.emp, rel.cov=rel.cov, f=2);
	llk.cov.rep1<-llk.cov.rep(llk1);

	llk.bal.rep1<-llk.bal.rep(d1, alpha.bal=2);
	llk.bal.err1<-llk.bal.err(d1, err=0.01);
	llk.bal.snp1<-llk.bal.snp(d1, err=0.01);

	op[i,1]<-i;
	op[i,2]<-llk1[1,2];

	op[i,3]<-op[i,2]+llk.bal.snp1;
	op[i,4]<-op[i,2]+llk.bal.err1;
	op[i,5]<-llk.cov.rep1+llk.bal.rep1;
	op[i,6]<-as.integer(d[i,2]=="REF_BUBBLE");

	a1<- d1[,1]>0;
	a2<- d1[,2]>0;
	max.gt<-rep(0, nrow(d1));
	max.gt[a1 & a2]<-1;
	max.gt[!a1 & a2]<-2;
	gt.cts<-hist(max.gt, plot=F, breaks=c(-0.5,0.5,1.5,2.5))$counts
	op[i,7]<-hwe(gt.cts, verb=F, ret=T);

	op[i,8]<-which.max(op[i,3:5]);
	if (llk1[1,2]<0) op[i,8]<-3;
	
}
cat("\n\nDone!\n\n");

thresh<-log(10);
min.cov<-8;
llk.thresh<- -40;
cv.cov.thresh<-10;
hwe.thresh<-0.05;
no.cov<-apply((d[1:to.do,2*ii+3]+d[1:to.do,2*ii+4])>0, 1, sum);
var.cov<-apply((d[1:to.do,2*ii+3]+d[1:to.do,2*ii+4]), 1, var);
cv.cov<-var.cov/apply((d[1:to.do,2*ii+3]+d[1:to.do,2*ii+4]), 1, mean);
#which.snp<-which((op[,8]==1) & ((op[,3]-apply(op[,4:5], 1, max))>thresh) & 
#	(no.cov>=min.cov) & op[,3]>llk.thresh & cv.cov<cv.cov.thresh & op[,7]>hwe.thresh);
which.snp<-which((op[,8]==1) & ((op[,3]-apply(op[,4:5], 1, max))>thresh) ) 


#cat("\nClassify ", 100*length(which.snp)/to.do, "% as SNPs", sep="");
#cat("\nOf those passing threshold, ", 100*mean(op[which.snp,6]), 
#	"% are REF BUBBLE compared to ", 100*mean(op[,6]), "% overall\n\n", sep="");
#which.fp<-which.snp[op[which.snp,6]==1];
#which.rep<-which(op[,8]==3);
#which.err<-which(op[,8]==2);
#cat("\n\nOf ", length(which.rep), " bubbles identified as repeat, ", 100*mean(op[which.rep,6]), 
#	"% are REF bubbles\n\n", sep="");
#cat("\n\nOf ", length(which.err), " bubbles identified as errors, ", 100*mean(op[which.err,6]), 
#	"% are REF bubbles\n\n", sep="");
#ss<-setdiff(which(op[,8]==1), which.snp);
#cat("\n\nOf ", length(ss), " bubbles identified as SNPs, but fail to meet threshold, ", 100*mean(op[ss,6]), 
#	"% are REF bubbles\n\n", sep="");
#cat("\n\nClassification eliminates ", 100-100*sum(op[which.snp,6])/sum(op[,6]), "% of errors\n\n", sep="");

#cts<-apply(d[,7:26], 1, sum);


## for each call, print variant name \t winning model \t confidence
for (i in 1:to.do) 
{
	cat(op[i,1], "\t", file="OUTPUT_FILENAME", append=TRUE,  sep=""); ##<<<<<<<<<<<< MODIFY - choose an output filename
	if (abs(op[i,8]-1)<0.5)
	{
		cat("variant\t", op[i,3]-max(op[i,4:5]),  file="OUTPUT_FILENAME", append=TRUE,  sep=""); ##<<<<<<<<<<<< MODIFY - choose an output filename
	}
	else if  (abs(op[i,8]-2)<0.5)
        {
               cat("error\t", op[i,4]-max(c(op[i,3],op[i,5])),	file="OUTPUT_FILENAME", append=TRUE,  sep=""); ##<<<<<<<<<<<< MODIFY - output filename
	}
	else if (abs(op[i,8]-3)<0.5)
	{
		cat("repeat\t", op[i,5]-max(op[i,3:4]),  file="OUTPUT_FILENAME", append=TRUE,  sep="");	 ##<<<<<<<<<<<< MODIFY - choose an output filename
	}
	else
	{
		cat("not classified", file="OUTPUT_FILENAME",append=TRUE, sep=""); ##<<<<<<<<<<<< MODIFY - choose an output filename
	}
	cat("\n", file="OUTPUT_FILENAME", append=TRUE, sep="");
}



