#===================================================================================
#  Lots of R packages for graphs
#	these are described in a CRAN task view here:
#	http://cran.r-project.org/web/views/gR.html
#	-install the 'network' package
#===================================================================================

#install.packages("network")
library(network)


#===================================================================================
#	Model neighborhood connections
#	-NEIGH: is a matrix of connections
#   -numneigh        #define the number of neighbors to connect
#   note: this leaves the ends with fewer connections
#===================================================================================
getNEIGH<-function(L,numneigh)
    {
    NEIGH<-matrix(0,L,L)  #set up an empty matrix
   

    for( i in numneigh:(L-numneigh))
        {
            for (o in 1:numneigh)
            {
        NEIGH[i,i-o]<-1
        NEIGH[i,i+o]<-1
            }
        }
        
    for(i in 1:numneigh)
        {
            for (o in 1:numneigh)
                {
                    NEIGH[i,o]<-1
                    NEIGH[(L+1)-i,(L+1)-o]<-1
                    NEIGH[(L+1)-i,o]<-1
                    NEIGH[i, (L+1)-o]<-1
                }
        }
    diag(NEIGH)<-0
    return(NEIGH)
    }

#NEIGH_NET<-network(NEIGH, directed=FALSE)

#===================================================================================
#	Model long-distance connections
#	-Each node has a number of connections drawn from  negative binomial distribution
#	-LD is the matrix of connections
#   -LD_NET is the network
#===================================================================================

getLD<-function(L,connectivity)
    {
        LD<-matrix(rbinom(L^2,1,connectivity),L)
            for( i in 1:L)
                {
                LD[i,]<-LD[,i]
                }
            diag(LD)<-0
        return(LD)
    }


#LD_NET<-network(LD, directed=FALSE)

#===================================================================================
#	create a single connection network out of NEIGH and LD
#   -ALL_NET is the matrix
#   -NET is the network
#===================================================================================
getALL_NET<-function(NEIGH,LD)
    {
    ALL_NET<-matrix(0,L,L)  #set up an empty matrix
    ALL_NET<-NEIGH+LD>0     #the unity of connections from long and neighborhood connections
    ALL_NET<-ALL_NET*1      #so that it is displayed in 0 and 1s
    return(ALL_NET)
    }

#NET<-network(ALL_NET, directed=FALSE)




#===================================================================================
#	Risk matrix, incorporating potential interventions
#       -RISK is the probability of infection among connection
#   Interventions:
#       -contain: an option to decrease risk for some columns which are assumed to report infection
#           -identification: % of cases identified in time to attempt containment
#           -participation: % of contacts of the node for whom risk will be reduced
#           -reduction: degree by which risk will be reduce, multiplicative
#       -cocoon (or any intervention trying to protect mother)
#           -cocoonable: % of mothers that could realistically be protected in time
#           -protection: Average reduction in risk to those protected (includes participation in protection efforts)
#===================================================================================


getrisk<-function(b, reduction=.4, identification=.5, participation=.6, cocoonable=.3, protection=.6, contain=FALSE, cocoon=FALSE)
    {
    RISK<-b*NEIGH
    
    if(contain==TRUE)
	   {
        identified<-sample(L, round(L*identification))    #randomly selects individuals that would present in time for ring treatment
        participants<-sample(L, round(L*participation))    #randomly selects those who would participate in ring treatment

            for (j in identified)
                {
                RISK[participants,j]<- RISK[participants,j]*reduction
                }
        }
    RISK<-RISK+(b*LD)       #add back the long distance risk that is not affected by containment
    
    if(cocoon==TRUE)
	   {
           cocoons<-sample(L, round(L*cocoonable))    #randomly selects individuals that could be protected
           
           for (j in cocoons)
           {
               RISK[,j]<- RISK[,j]*protection
           }
       }
    
    return(RISK)
    }



#===================================================================================
#	Plotting functions
#	-COORDS lets it draw network once and saves the coordinates of each node
#	-COLS a vector to define colors in the plot 
#	-visualize: plots network, coloring the I nodes red
#===================================================================================




visualize_net<-function(NET)
	{
    COORDS<<-plot(NET,vertex.col="black")
    COLS<<-rep(1,L)
	COLS[I==1]<-"red"
	COLS[E==1]<-"orange"
	COLS[S==1]<-"black"
	COLS[R==1]<-"white"
	plot(NET,vertex.col=COLS, jitter=FALSE, coord=COORDS)
	}


visualize_nodes<-function()
	{
	COLS[I==1]<-"red"
	COLS[E==1]<-"orange"
	COLS[S==1]<-"black"
	COLS[R==1]<-"white"
    lines(COORDS,col=COLS,pch=19,type="p") # much faster than points()
	}

#===================================================================================
#	-Set up Simulation 
#	-S,E,I,R vectors
#===================================================================================

setup<-function()
	{
	S<<-rep(1,L)
	E<<-rep(0,L)
	I<<-rep(0,L)
	R<<-rep(0,L)
	recoverday<<-rep(0,L)
    cases<<-NA
	}

#===================================================================================
#	Key functions of stochastic simulation
#	-infect changes I from 0 to 1, also changes E and S to 0
#===================================================================================

infect<-function(node, day=i)
	{
        I[node]<<-1
        E[node]<<-0
        S[node]<<-0
        cases<<-which(I==1)
        recoverday[node]<<-day+duration
	}
		
expose<-function()
	{
        if(length(cases==1))
           {
               temp<-ALL_NET[,cases]==1
           }
           
        if(length(cases)>1)
           {
               temp<-rowSums(ALL_NET[,cases])>=1
           }
        if(length(temp>=1))
            {
                temp <- temp & !R # recovered are not exposed assuming immunity at least over time period
                E[temp]<<-1
                E[cases]<<-0
                S[temp]<<-0
            }
	}

recover<-function(node)
	{
        I[node]<<-0
        E[node]<<-0
        S[node]<<-0
        R[node]<<-1
        cases<<-which(I==1)
	}


#===================================================================================
#	Simulation function
#   Parameters:
#       -b : probaility infection given exposure / time step (ie rate)
#       -endtime: number of time steps in simulation
#       -index_cases: the # of cases at timepoint 1
#       -connectivity % of population that comes into relevant contact. Note this scales with L, the population size
#       -numneigh Number of neighboring nodes a cases will expose, also the number that might be treated in containment. Currently the number is doubled as the simulation is 2 dimensional (a 3 means 6 nodes will be exposed)
#       -coco : whether a cocoon strategy is attempted
#       -cont : whether a containment strategy is attempted
#       -identification, participation, reduction, coconnable, protection (see the Risk Matrix function above)
#===================================================================================

dosim<-function(cont= FALSE, coco=FALSE, index_cases=1, numneigh=3, connectivity=0.05, b=.01, identification=1, participation=1, reduction=0,cocoonable=.1, protection=.6)
    {
        #reset the population vectors
            i<<-1
            setup()
            PREV<<-temp<<-rep(0,endtime)  #empty vector
        
        #create the risk matrices
            RISK<<-getrisk(b, reduction=reduction, identification=identification, participation=participation, cocoonable=cocoonable, protection=protection, contain=cont, cocoon=coco)
            
        #seed the epidemic
        
            cases<<-(sample(L,index_cases))           # infect the index cases
            infect(cases)
            expose()                                # expose nodes connected to index cases
       
        #simulate the epidemic
           for(i in 2:endtime)
                {
                    
                    random<-matrix(runif(L*L),L,L)  # draw a random number matrix for stochastic infection
                    risk<-E*RISK                    # multiplying risk matrix by the exposed vector
                    toinfect<-which(rowSums(random<risk)>0)
                    infect(toinfect,day=i)
                    recover(which(recoverday==i))
                    expose()
                    PREV[i]<-sum(I)/L               # store the prevalence at time i
                }
        #return total cases and prevalence over time
            totcases<-sum(R)+sum(I)
            togo<-c(PREV, totcases)
            return(togo)
    }


#===================================================================================
#   Global Parameters
#       -sims: number of simulations per arm
#       -duration : number of timepoints a case is infectious
#
#===================================================================================
sims<-10
endtime<-100
L<-200
I<-rep(0,L)     #allows use of the symbol I for infectious and sets it to zero
duration<-1
numneigh<-3
connectivity<-0.05

baseline<-matrix(NA,endtime,sims)
contained<-matrix(NA,endtime,sims)
cocooned<-matrix(NA,endtime,sims)
both<-matrix(NA,endtime,sims)

basecase<-rep(NA,sims)
contcase<-rep(NA,sims)
cococase<-rep(NA,sims)
bothcase<-rep(NA,sims)


for (x in 1:sims)
    {
    #create the exposure and risk matrices
        NEIGH<<-getNEIGH(L, numneigh)
        LD<<-getLD(L, connectivity)
        ALL_NET<<-getALL_NET(NEIGH,LD)
        

    #no containment
        temp<-dosim(cont= FALSE, coco=FALSE)
        baseline[,x]<-temp[1:endtime]
        basecase[x]<-temp[endtime+1]
    #containment

        temp<-dosim(cont= TRUE, coco=FALSE)
        contained[,x]<-temp[1:endtime]
        contcase[x]<-temp[endtime+1]

    #cocoon

        temp<-dosim(cont= FALSE, coco=TRUE)
        cocooned[,x]<-temp[1:endtime]
        cococase[x]<-temp[endtime+1]
    #both

        temp<-dosim(cont= TRUE, coco=TRUE)
        both[,x]<-temp[1:endtime]
        bothcase[x]<-temp[endtime+1]
    }

#===================================================================================
#	Visualize the results
#	-black lines: simulations with no interventions
#	-blue lines: Simulaitons with containment
#   -purple lines: simulations with cocooning
#   -red lines: simulations with both containment and cocooning
#
#   -Thicker lines are means across the simulations in each arm
#===================================================================================


#plots of each simulation
quartz()
matplot(1:endtime, baseline, type = "l", xlab = "Time", ylab = "Prevalence", main = "", lwd = .1, lty = 1, bty = "l", col = "black")
matplot(1:endtime, contained, type = "l", xlab = "Time", ylab = "Prevalence", main = "", lwd = .1, lty = 1, bty = "l", col = "blue", add=TRUE)
matplot(1:endtime, cocooned, type = "l", xlab = "Time", ylab = "Prevalence", main = "", lwd = .1, lty = 1, bty = "l", col = "purple", add=TRUE)
matplot(1:endtime, both, type = "l", xlab = "Time", ylab = "Prevalence", main = "", lwd = .1, lty = 1, bty = "l", col = "red", add=TRUE)

#show the means for each set of simulations
lines(rowMeans(baseline),col="black", lwd=3)
lines(rowMeans(contained),col="blue", lwd=3)
lines(rowMeans(cocooned),col="purple", lwd=3)
lines(rowMeans(both),col="red", lwd=3)
legend(60, max(baseline), c("baseline", "containment", "cocooning", "both"), pch = 18, col = c("black", "blue", "purple", "red"))

#show final epidemic size
basecase
contcase
cococase
bothcase




#===================================================================================
#	-par sets up a plot window. par(ask=TRUE) requires a <Enter> between plots (for a movie)
#   -Visualize_net shows the network of local and longdistnace contacts
#===================================================================================
dosim(connectivity=0)
par(ask=FALSE)
visualize_net(network(ALL_NET, directed=FALSE))



#####-------------old below

#===================================================================================
#	-Simulated population parameters
#===================================================================================

L<-num_obs<-100  #number of households and individuals
S<-rep(1,L)
E<-rep(0,L)
I<-rep(0,L)
R<-rep(0,L)





#define the length in days of the simulation
reps<-100

#set time to day 1
i<-1

#an empty vector to store the prevalence over the course of the simulation
PREV<-rep(0,reps)

#sets the whole populations to susceptible
setup()

#set how many index cases and draw them randomly
index<-index_case<-sample(L,5)

#infect the index cases
infect(index_case)

#expose those nodes connected to the index cases
expose()



