/*********************************************
 * OPL 12.8.0.0 Model
 * Author: sosa01
 * Creation Date: 16.11.2020 at 17:04:24
 *********************************************/
using CP;

int timeLimit = ...;

int numberWorkers = ...;
int numberProducts = ...;
int numberVariants = ...;
int numberStations = ...;
range allWorkers = 1..numberWorkers;
range allProducts = 1..numberProducts;
range allVariants = 1..numberVariants;
range allStations = 1..numberStations;

tuple workerTuples {key int id; int skill; int releasetime; int currentStation;}; 
{workerTuples} workers = ...;

tuple productTuples {key int id; int variant; int releasetime; int currentStation; int currentBuffer; int currentStep;};
{productTuples} products = ...;

int workerWalkingTimes[allStations][allStations] = ...;
tuple triplet {int id1; int id2; int value;}; 
{triplet} workerTransitionTimes;

int workstepDurations[allVariants][allStations] = ...;
int isAutomatic[allVariants][allStations] = ...;
int setupTimes[allVariants][allStations] = ...;
int takedownTimes[allVariants][allStations] = ...;

//dvar int productionTimeStation[allProducts][allStations];
//dvar int productionTime[allProducts];
//dvar int firstStationTime[1..numberProducts-1];


dvar interval work[worker in allWorkers][product in 0..numberProducts][station in allStations][workstep in 1..2] optional;
dvar interval pro[product in allProducts][station in allStations][workstep in 1..2] optional;

dvar sequence worker_seq[worker in allWorkers] in all(product in 0..numberProducts, station in allStations, workstep in 1..2) work[worker][product][station][workstep]
						types all(product in 0..numberProducts, station in allStations, workstep in 1..2) station;

// set time limit
execute {
	cp.param.timeLimit = timeLimit;
	for(var i in allStations){
		for(var j in allStations){
			workerTransitionTimes.add(i, j, workerWalkingTimes[i][j]);	
		}	
	}
}

minimize max(k in 1..2) endOf(pro[numberProducts][numberStations][k]);

subject to{

//production time
forall(worker in workers)
  forall(product in products){ 
  // product 0 is a dummy product for the initial position of the worker
  	(presenceOf(work[worker.id][0][worker.currentStation][1])) =>
    	sizeOf(work[worker.id][0][worker.currentStation][1]) == 0;
    forall(station in allStations){
    	if(isAutomatic[product.variant][station] == 1){
    		(presenceOf(work[worker.id][product.id][station][1])) =>
    			sizeOf(work[worker.id][product.id][station][1]) == round((100/worker.skill) * setupTimes[product.variant][station]);
    		(presenceOf(work[worker.id][product.id][station][2])) =>
    			sizeOf(work[worker.id][product.id][station][2]) == round((100/worker.skill) * takedownTimes[product.variant][station]);
    	} 
    	else{
    		(presenceOf(work[worker.id][product.id][station][1])) =>
    		sizeOf(work[worker.id][product.id][station][1]) == round((100/worker.skill) * workstepDurations[product.variant][station]);
    	}
  	}
  }  

  
forall(product in products){
	// for each product all production steps (from the initial station onwards) at a station which have time > 0 must be present
	// for the initial station only production steps that have not occurred yet must be produced
	if(product.currentStation > 0){
		//all steps before the initial step are not present	
		forall(station in 1..product.currentStation-1, workstep in 1..2){
			presenceOf(pro[product.id][station][workstep]) == 0;
		}
		//all steps after the initial station are present
		forall(station in product.currentStation+1..numberStations){
			presenceOf(pro[product.id][station][1]) == 1;
			presenceOf(pro[product.id][station][2]) == isAutomatic[product.variant][station];		 			
		}
		presenceOf(pro[product.id][product.currentStation][1]) == 0;
		//at the initial station only steps after the initial step are present
		if(isAutomatic[product.variant][product.currentStation] == 1 && product.currentStep != 4){
			presenceOf(pro[product.id][product.currentStation][2]) == 1;
		}
		else{
			presenceOf(pro[product.id][product.currentStation][2]) == 0;		
		}
	}
	//if no workstep has been interrupted and the product is either on a buffer or the initial queue	
	else{
		forall(station in 1..product.currentBuffer, workstep in 1..2){
			presenceOf(pro[product.id][station][workstep]) == 0;	
		}
		forall(station in product.currentBuffer+1..numberStations){
		 	presenceOf(pro[product.id][station][1]) == 1;
		 	presenceOf(pro[product.id][station][2]) == isAutomatic[product.variant][station];
		}
	}
}


//each production Step only occurs once
forall(product in products)
  forall(station in allStations)
    forall(workstep in 1..2){
    	alternative(pro[product.id][station][workstep], all(worker in allWorkers) work[worker][product.id][station][workstep]);
    	synchronize(pro[product.id][station][workstep], all(worker in allWorkers) work[worker][product.id][station][workstep]);
}


//precedence for stations
forall(product in products){
  forall(station1,station2 in allStations: station1<station2)
    forall(workstep1,workstep2 in 1..2){
   		endBeforeStart(pro[product.id][station1][workstep1], pro[product.id][station2][workstep2]);
  }
}

//precedence for products
forall(product1,product2 in products: product1.id<product2.id)
  forall(station in allStations)
    forall(workstep1,workstep2 in 1..2){
    	endBeforeStart(pro[product1.id][station][workstep1], pro[product2.id][station][workstep2]);
  }   

//buffer
forall(product1,product2 in products: product1.id<product2.id){
	forall(station in 1..numberStations-1){
		 startBeforeEnd(pro[product1.id][station+1][1], pro[product2.id][station][isAutomatic[product2.variant][station]+1]);
	}    
}	

//walking time/ precedence for workers
forall(workstep in allWorkers){
		noOverlap(worker_seq[workstep], workerTransitionTimes, true);
		noOverlap(worker_seq[workstep]);
}

//automatic step time
forall(product in products)
  	forall(station in allStations){
		(isAutomatic[product.variant][station] == 1) 
			=> endOf(pro[product.id][station][1]) + workstepDurations[product.variant][station] == startOf(pro[product.id][station][2]);
}


// initial setup
forall(worker in workers){
	// if the worker is still finishing a job, then all other jobs he does must start after the releasetime (plus the walking time to the other station)
	forall(product in 0..numberProducts)
		forall(station in allStations)
		    forall(workstep in 1..2){
				if (worker.releasetime > 0) {
					presenceOf(work[worker.id][product][station][workstep])
						=> (startOf(work[worker.id][product][station][workstep]) >= worker.releasetime + workerWalkingTimes[worker.currentStation][station]);		      
		      	}
	}
	// the worker must start working on product "0" for time 0 at the workers initial station
	presenceOf(work[worker.id][0][worker.currentStation][1]) == 1;
	first(worker_seq[worker.id], work[worker.id][0][worker.currentStation][1]);
	// and no other intervals of product 0 exist for this worker	
	forall(station in allStations)
	  forall(workstep in 1..2){
	    (station != worker.currentStation || workstep != 1) => presenceOf(work[worker.id][0][station][workstep]) == 0;
	}
}



// if a product is still being worked on at a station, this station is "blocked" until the release time
forall(product, product2 in products: product.id < product2.id){
  	if(product.releasetime > 0){
  		if(product.currentStation != 0 && (isAutomatic[product.variant][product.currentStation] == 0 || product.currentStep == 4)){
  			presenceOf(pro[product2.id, product.currentStation, 1]) => startOf(pro[product2.id, product.currentStation, 1]) >=  product.releasetime;		
  		}
  	}  
}

forall(product in products){
	if(product.releasetime > 0 && product.currentStation != 0 && isAutomatic[product.variant][product.currentStation] == 1){
		if(product.currentStep == 2){
			presenceOf(pro[product.id, product.currentStation, 2]) => startOf(pro[product.id, product.currentStation, 2]) >= product.releasetime + workstepDurations[product.variant][product.currentStation];
		}
		else if(product.currentStep == 3){
			presenceOf(pro[product.id, product.currentStation, 2]) => startOf(pro[product.id, product.currentStation, 2]) >= product.releasetime;
		}
	}	
}


//additional duplicate constraints that should help speed up the search

forall(worker in allWorkers){
	last(worker_seq[worker], work[worker][numberProducts][numberStations][2]);
}
    
}

tuple solution
{
  int worker;
  int product;
  int station;
  int workstep;
  int start;
  int end;
  int exists;
}

{solution} solutions={<w, p, s, k, startOf(work[w][p][s][k]), endOf(work[w][p][s][k]), presenceOf(work[w][p][s][k])> | w in allWorkers,
				 p in allProducts, s in allStations, k in 1..2};

 