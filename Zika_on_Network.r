#===================================================================================
#  Lots of R packages for graphs
#	these are described in a CRAN task view here:
#	http://cran.r-project.org/web/views/gR.html
#	-install the 'network' package
#===================================================================================

#install.packages("network")
library(network)
#help("network-package")


#===================================================================================
#	-Simulated population parameters
#===================================================================================

L<-num_obs<-100  #number of households and individuals
S<-rep(1,L)
E<-rep(0,L)
I<-rep(0,L)
R<-rep(0,L)


#===================================================================================
#	Model neighborhood connections
#	-NEIGH: is a matrix of connections
#	-NEIGH_NET is the same saved as a network
#   note: this leaves the ends with fewer connections
#===================================================================================

NEIGH<-matrix(0,L,L)  #set up an empty matrix
numneigh<-3        #define the number of neighbors to connect

for( i in numneigh:L-numneigh)
	{
        for (o in 1:numneigh)
        {
    NEIGH[i,i-o]<-1
    NEIGH[i,i+o]<-1
        }
    }
diag(NEIGH)<-0

NEIGH_NET<-network(NEIGH, directed=FALSE)

#===================================================================================
#	Model long-distance connections
#	-Each node has a number of connections drawn from  negative binomial distribution
#	-LD is the matrix of connections
#   -LD_NET is the network
#===================================================================================

connectivity<-.01

LD<-matrix(rbinom(L^2,1,connectivity),L)
for( i in 1:L)
    {
    LD[i,]<-LD[,i]
    }

diag(LD)<-0

LD_NET<-network(LD, directed=FALSE)

#===================================================================================
#	create a single connection network out of NEIGH and LD
#   -ALL_NET is the matrix
#   -NET is the network
#===================================================================================

ALL_NET<-NEIGH+LD>0  #the unity of connections from long and neighborhood connections
ALL_NET<-ALL_NET*1       #so that NET is displayed in 0 and 1s
NET<-network(ALL_NET, directed=FALSE)




#===================================================================================
#	Risk matrix, incorporating potential interventions
#       -RISK is the probability of infection among connection
#   Interventions:
#       -contain: an option to decrease risk for some columns which are assumed to report infection
#           -identification: % of cases identified in time to attempt containment
#           -participation: % of contacts of the node for whom risk will be reduced
#           -reduction: degree by which risk will be reduce, multiplicative
#===================================================================================



getrisk<-function(b, contain=FALSE, cocoon=FALSE)
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
    return(RISK)
    }



#===================================================================================
#	Plotting functions
#	-COORDS lets it draw network once and saves the coordinates of each node
#	-COLS a vector to define colors in the plot 
#	-visualize: plots network, coloring the I nodes red
#===================================================================================


COORDS<-plot(NET,vertex.col="black")
COLS<-rep(1,L)

visualize_net<-function()
	{
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
	}

#===================================================================================
#	Key functions of stochastic simulation
#	-infect changes I from 0 to 1, also changes E and S to 0
#===================================================================================

infect<-function(node,day=i)
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
	   index<-ALL_NET[,cases]==1
	   }
       
	if(length(cases)>1)
	   {
	   index<-rowSums(ALL_NET[,cases])>=1
	   }
       
	index<- index & !R # recovered are not exposed assuming immunity at least over time period
	E[index]<<-1
	E[cases]<<-0
	S[index]<<-0
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
#	Simulation setup
#	-reps: define how many time steps
#	-setup()  set S to 1, I, E, R to 0
#	-PREV: an empty vector to keep track of prevalence at each timestep
#	-par sets up a plot window. par(ask=TRUE) requires a <Enter> between plots (for a movie)
#	-b : probaility infection given exposure / time step (ie rate)
#	-duration : number of days infectious
#	-Assign index cases
#	-infect it with infect()
#	-expose its neighbors and long distance contacts with expose()
#===================================================================================

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

#plot the starting conditions
par(ask=FALSE)
visualize_net()



#probability of infection in exposed node per time step
b <-.1

#duration of infectiousness
duration<-5

#containment parameters
identification<-.5
participation<-.6
reduction<-.4


#create the risk matrix
RISK<-getrisk(b, contain=TRUE, cocoon=FALSE)




#===================================================================================
#	Watch a single simulation
#	-reps: define how many time steps
#	-outputs a movie 
#	-plots prevalence over time
#===================================================================================
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


#probability of infection in exposed node per time step
b <-.1

#duration of infectiousness
duration<-5

#containment parameters
identification<-.5
participation<-.6
reduction<-.4


#create the risk matrix
RISK<-getrisk(b, contain=TRUE, cocoon=FALSE)


#plot the starting conditions
par(ask=FALSE)
visualize_net()

    for(i in 2:reps)
    {
        random<-matrix(runif(L*L),L,L)  #draw a random number matrix for stochastic infection
        risk<-E*RISK		#multiplying risk matrix by the exposed vector
        toinfect<-which(rowSums(random<risk)>0)
        infect(toinfect)
        recover(which(recoverday==i))
        expose()
        visualize_nodes()   #plot the current state of the network
        PREV[i]<-sum(I)/L   #store the prevalence at time i
        #CUMULPREV <- sum(I+R)/L
    }

plot(PREV,typ="l", ylab="Prevalence",xlab="Time step",ylim=c(0,1))




#===================================================================================
#	Simulation loop
#	-reps: define how many time steps
#	-outputs a movie (currently commented out)
#	-plots prevalence over time
#===================================================================================

dosim<-function()
{
    for(i in 2:reps)
    {
        random<-matrix(runif(L*L),L,L)  #draw a random number matrix for stochastic infection
        risk<-E*RISK		#multiplying risk matrix by the exposed vector
        toinfect<-which(rowSums(random<risk)>0)
        infect(toinfect)
        recover(which(recoverday==i))
        expose()
        #visualize_nodes()   #plot the current state of the network
        PREV[i]<-sum(I)/L   #store the prevalence at time i
        return(PREV)
        #CUMULPREV <- sum(I+R)/L
    }
}

sims<-60
Ps<-matrix(NA,reps,sims)
for (x in 1:sims)
{
    Ps[,x]<-dosim()
}

matplot(1:reps, Ps, type = "l", xlab = "Time", ylab = "Prevalence", main = "", lwd = 1, lty = 1, bty = "l", col = "blue")


plot(PREV,typ="l", ylab="Prevalence",xlab="Time step",ylim=c(0,1))









#---------old below


#===================================================================================
#	Simple Simulation loop
#	-reps: define how many time steps
#	-outputs a movie (currently commented out)
#	-plots prevalence over time
#===================================================================================

for(i in 2:reps)
{
    random<-runif(L)
    risk<-E*b		#multiplying the risk by whether or not exposed
    infect(random<risk)
    recover(which(recoverday==i))
    expose()
    visualize_nodes()   #plot the current state of the network
    PREV[i]<-sum(I)/L   #store the prevalence at time i
    #CUMULPREV <- sum(I+R)/L
}

plot(PREV,typ="l", ylab="Prevalence",xlab="Time step",ylim=c(0,1))



