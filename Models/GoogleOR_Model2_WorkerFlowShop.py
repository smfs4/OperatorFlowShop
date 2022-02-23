class CircuitORSolver(object):
    """Solver Class of the JobShop Problem"""

    def CircuitBoschJobshopSat(self, config_path, init_path, sol_path, max_time):
        """Minimal jobshop problem."""
        import collections
        import time
        import xml.etree.ElementTree as ET
        import numpy as np
        import struct

        from XML_Reader import XML_Reader
        from XML_Writer import XML_Writer


        # Import Python wrapper for or-tools CP-SAT solver.
        from ortools.sat.python import cp_model

        problem = XML_Reader(config_path, init_path)
        production_time_matrix = problem.get_production_time_matrix()
        worker_movement_time_matrix = problem.get_worker_movement_matrix()
        automatic_worksteps_times_matrix = problem.get_automatic_times_matrix()
        number_manual_worksteps = problem.get_number_worksteps()
        releasetime_workers = problem.get_releasetime_workers()
        leftover_job_products = problem.get_leftover_jobs_products()
        worker_skill = problem.get_worker_skill()
        products_type = problem.get_products_type()
        products_initial = problem.get_products_initial()
        workers_initial_station = problem.get_worker_initial_stat()


        tic = time.perf_counter()
        model = cp_model.CpModel()

        stations_count = len(worker_movement_time_matrix[0])
        print("Number of stations is", stations_count)
        all_stations = range(stations_count)

        products_count = len(products_type)
        print("Number of products is", products_count)
        all_products = range(products_count)

        num_product_types = 4 #len(production_time_matrix)
        print("Number of product types is", num_product_types)

        horizon = 1000000
        print("The horizon is", horizon)

        workers_count = len(worker_skill)
        print("The number of workers is", workers_count)
        all_workers = range(workers_count)

        # Named tuple to store information about created variables.
        #task_type = collections.namedtuple('task_type', 'start end duration interval')

        # Named tuple to store information about created variables.
        opt_task_type = collections.namedtuple('task_type', 'exists start end duration interval')
        task_type = collections.namedtuple('task_type', 'start end duration interval')

        # Named tuple to manipulate solution information.
        assigned_task_type = collections.namedtuple('assigned_task_type',
                                                    'start exists product station workstep duration')

        # Creates job intervals, for two intervals assigned to the same worker, whether one occurs before
        # the other(is_before), whether both intervals exist for the same worker (both_exist) and the start_time_diff
        # (the difference in start times between two jobs if they both exist for this worker, otherwise we set this to 0)
        # and add to the corresponding machine lists.
        all_tasks = {}
        all_tasks_workers = {}
        node_name = {}
        arc_workers = [[] for worker in all_workers] 
        arc_literals = {}
        is_first = {}
        task_starts = [] ################################################

        counter = 1
        for product in all_products:
            for station in range(products_initial[0][product], stations_count):
                for workstep in range(number_manual_worksteps[products_type[product]][station]):
                    if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                        continue
                    suffix = '_%i_%i_%i' % (product, station, workstep)
                    #duration of each job
                    duration_var = model.NewIntVar(0, horizon, 'duration' + suffix)
                    #start time of each job
                    start_var = model.NewIntVar(0, horizon, 'start' + suffix)
                    #end time of each job
                    end_var = model.NewIntVar(0, horizon, 'end' + suffix)
                    #interval variable with these start, end, duration variables
                    interval_var = model.NewIntervalVar(start_var, duration_var, end_var,
                                                        'interval' + suffix)
                    #dictionary of all tasks
                    all_tasks[product, station, workstep] = task_type(
                        start=start_var, end=end_var, duration = duration_var, interval=interval_var)
                    node_name[product, station, workstep] = counter
                    counter+=1
                    task_starts.append(start_var) #############################################
                    for worker in all_workers: 
                        is_first_var = model.NewBoolVar('is_first' + suffix)
                        is_last_var = model.NewBoolVar('is_last' + suffix)
                        arc_workers[worker].append((0, node_name[product, station, workstep], is_first_var))
                        arc_workers[worker].append((node_name[product, station, workstep], 0, is_last_var))
                        is_first[worker, product, station, workstep] = is_first_var
                        # for each worker we create an optional interval for each job
                        suffix_worker = '_%i_%i_%i_%i' % (worker, product, station, workstep)
                        # variable saying whether the given job is produced by the worker
                        exists_var = model.NewBoolVar('exists' + suffix_worker)
                        # duration for the jobs
                        duration_worker = int(round(float(production_time_matrix[products_type[product]][station][workstep]) * (float(100)/float(worker_skill[worker]))))
                        # interval variable with the same start and end time as the "general" job and the duration specific
                        # to the worker (ok since only ever one of these interval variables exist for each job)
                        # dictionary of optional intervals
                        all_tasks_workers[product, station, workstep, worker] = opt_task_type(
                            exists = exists_var, start=start_var, end=end_var, duration = duration_worker, interval=interval_var)
                        #creates a defaultdict of all interval variables for each worker
                        #worker_to_intervals[worker].append(interval_var_worker)
                        arc_literals[product, station, workstep, product, station, workstep, worker] = exists_var.Not()
                        arc_workers[worker].append([node_name[product, station, workstep], node_name[product, station, workstep], exists_var.Not()])

        print("Number of operations is", counter)
                                    
        # define the arcs between the workstep nodes for each worker
        for product in all_products:
            for station in range(products_initial[0][product], stations_count):
                for workstep in range(number_manual_worksteps[products_type[product]][station]):
                    if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                        continue
                    for worker in all_workers:    
                        for product2 in all_products:
                            for station2 in range(products_initial[0][product2], stations_count):
                                for workstep2 in range(number_manual_worksteps[products_type[product2]][station2]):
                                    if (station == station2) and (product == product2) and (workstep == workstep2):
                                        continue
                                    suffix_is_before = '_%i_%i_%i_%i_%i_%i_%i' % (worker, product, station, workstep, product2, station2, workstep2)
                                    # creates a boolean variable which is true if worker produces both (product, station) and (product2, station2)
                                    arc_literal_var = model.NewBoolVar('arc_literal' + suffix_is_before)
                                    arc_literals[product, station, workstep, product2, station2, workstep2, worker] = arc_literal_var
                                    arc_workers[worker].append([node_name[product, station, workstep], node_name[product2, station2, workstep2], arc_literal_var])
  
        # each interval is present for exactly one worker
        print("1")
        for product in all_products:
            for station in range(products_initial[0][product], stations_count):
                for workstep in range(number_manual_worksteps[products_type[product]][station]):
                    if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                        continue
                    model.Add(sum((all_tasks_workers[product, station, workstep, worker].exists)
                                            for worker in all_workers) == 1)

        # if the release time for a worker is greater than 0, then he can't start working until then
        print("2")
        for worker in all_workers:
            if (releasetime_workers[worker] > 0):
               for product in all_products:
                    for station in range(products_initial[0][product], stations_count):
                        for workstep in range(number_manual_worksteps[products_type[product]][station]):
                            if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                                continue
                            model.Add(all_tasks_workers[product, station, workstep, worker].start >= releasetime_workers[worker] +
                                     worker_movement_time_matrix[workers_initial_station[worker]][station]).OnlyEnforceIf(
                                    is_first[worker, product, station, workstep])
        
        # if a product has some leftover time at a station, the product cannot start being processed elsewhere before the releasetime
        # and no other product can be produced at that station until the product is done
        print("3")
        for product in all_products:
            for station in range(products_initial[0][product], stations_count):
                if (leftover_job_products[0][product] <= 0):
                    continue
                for workstep in range(number_manual_worksteps[products_type[product]][station]):
                    if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                        continue
                    model.Add(all_tasks[product, station, workstep].start >= leftover_job_products[0][product])
                # the station is blocked for all further products until the release time
                for product2 in range(product+1, products_count):
                    model.Add(all_tasks[product2, leftover_job_products[1][product]-1, 0].start >= leftover_job_products[0][product])
                # the next buffer must be clear when the product finishes (but only if the final workstep at the station is being completed before the release time)
                if products_initial[1] == 1:
                    for product1 in range(0, product):
                        if products_initial[0][product1] == leftover_job_products[1][product]:
                            continue
                        model.Add(all_tasks[product1, leftover_job_products[1][product], 0].start <= leftover_job_products[0][product])
                        
        #precedence constraints stations
        print("4")
        for product in all_products:
            for station1 in range(products_initial[0][product], stations_count):
                for station2 in range(station1+1,stations_count):
                    model.Add(all_tasks[product, station1, number_manual_worksteps[products_type[product]][station1]-1].end <= all_tasks[product, station2, 0].start)

        #precedence constraints products
        print("5")
        for product1 in all_products:
            for product2 in range(product1+1, products_count):
                for station in range(products_initial[0][product1], stations_count):
                    if (station == products_initial[0][product2] and products_initial[1][product2] != 1):
                        continue
                    model.Add(all_tasks[product1, station, number_manual_worksteps[products_type[product1]][station]-1].end <= all_tasks[product2, station, 0].start)


        #buffer constraint
        print("6")
        for product1 in all_products:
            for product2 in range(product1+1, products_count):
                for station in range(max([products_initial[0][product1]-1, products_initial[0][product2]]), stations_count-1):
                    if (station+1 == products_initial[0][product1] and products_initial[1][product1] != 1):
                        continue
                    model.Add(all_tasks[product1, station+1, 0].start <= all_tasks[product2, station, number_manual_worksteps[products_type[product2]][station]-1].end)
  
    
        #automatic station
        print("7")
        for product in all_products:
            for station in range(products_initial[0][product], stations_count):
                if (station == products_initial[0][product] and 0 < products_initial[1][product]-1):
                    continue
                if (number_manual_worksteps[products_type[product]][station] == 2):
                    model.Add(all_tasks[product, station, 0].end == all_tasks[product, station, 1].start - automatic_worksteps_times_matrix[products_type[product]][station])
        
        # circuit and walking constraints
        print("8")
        for worker in all_workers:
            model.AddCircuit(arc_workers[worker])
            for product1 in all_products:
                for product2 in all_products:
                    for station1 in range(products_initial[0][product1], stations_count):
                        for station2 in range(products_initial[0][product2], stations_count):
                            for workstep1 in range(number_manual_worksteps[products_type[product1]][station1]):
                                if station1 == products_initial[0][product1] and (workstep1 < products_initial[1][product1]-1):
                                    continue
                                for workstep2 in range(number_manual_worksteps[products_type[product2]][station2]):
                                    if station2 == products_initial[0][product2] and (workstep < products_initial[1][product2]-1):
                                        continue
                                    if (station1 != station2) or (product1 != product2) or (workstep1 != workstep2):
                                        model.Add(all_tasks[product2, station2, workstep2].start >= 
                                                        (all_tasks[product1, station1, workstep1].end + 
                                                        worker_movement_time_matrix[station2][station1])).OnlyEnforceIf(
                                              arc_literals[product1, station1, workstep1, product2, station2, workstep2, worker])
                                    

        # Initial positions of the workers means that no production step can start before the worker has
        # walked to the station from its initial position
        print("9")
        for worker in all_workers:
            for product in all_products:
                for station in range(products_initial[0][product], stations_count):
                    for workstep in range(number_manual_worksteps[products_type[product]][station]):
                        if station == products_initial[0][product] and (workstep < products_initial[1][product]-1):
                            continue
                        model.Add(all_tasks_workers[product, station, workstep, worker].start >= 
                                worker_movement_time_matrix[workers_initial_station[worker]][station]).OnlyEnforceIf(
                                    all_tasks_workers[product, station, workstep, worker].exists)


                
        # Makespan objective.
        obj_var = model.NewIntVar(0, horizon, 'makespan')
        model.Add(obj_var == all_tasks[products_count-1, stations_count-1, number_manual_worksteps[products_type[products_count-1]][stations_count-1]-1].end)
        model.Minimize(obj_var)

        toc = time.perf_counter()
        print(f"Solver called after {toc - tic:0.4f} seconds")
        tada = time.perf_counter()
        model_loading_time = tada - tic
      
        # Solve model.
        model_loading_time = tada - tic
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 16# if multiple workers
        solver.parameters.max_time_in_seconds = max_time - model_loading_time
        status = solver.Solve(model)
        toc = time.perf_counter()
        print("solved")
        if status == cp_model.INFEASIBLE:
            print(f"Shown infeasible in {toc - tic:0.4f} seconds")
            status_message = "INFEASIBLE"

        if status == cp_model.MODEL_INVALID:
            print(f"Shown model invalid in {toc - tic:0.4f} seconds")
            status_message = "MODEL_INVALID"

        if status == cp_model.OPTIMAL:
            print(f"Solved optimally in {toc - tic:0.4f} seconds")
            status_message = "OPTIMAL"
            # Create one list of assigned tasks per machine.
        if status == cp_model.FEASIBLE:
            print(f"Solved suboptimally in {toc - tic:0.4f} seconds")
            status_message = "SUBOPTIMAL"
      
        if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
            # Finally print the solution found.
            solution_objective = solver.ObjectiveValue()
            solution_time = toc - tic
            solution_status = status_message
            print('Schedule Length: %i' % solution_objective)

            assigned_stations = collections.defaultdict(list)
            for product_id in all_products:
                for station_id in all_stations:
                    if station_id < products_initial[0][product_id]:
                        continue
                    for workstep_id in range(number_manual_worksteps[products_type[product_id]][station_id]):
                        if station_id == products_initial[0][product_id] and (workstep_id < products_initial[1][product_id]-1):
                            continue
                        assigned_stations[product_id].append(
                            assigned_task_type(
                                start=solver.Value(all_tasks[product_id, station_id, workstep_id].start),
                                exists = True,
                                product=product_id,
                                station=station_id,
                                workstep=workstep_id,
                                duration=solver.Value(all_tasks[product_id, station_id, workstep_id].duration)))

            # Create one list of assigned tasks per worker
            assigned_products_workers = collections.defaultdict(list)
            for worker_id in all_workers:
                for product_id in all_products:
                    for station_id in all_stations:
                        if station_id < products_initial[0][product_id]:
                            continue
                        for workstep_id in range(number_manual_worksteps[products_type[product_id]][station_id]):
                            if station_id == products_initial[0][product_id] and (workstep_id < products_initial[1][product_id]-1):
                                continue
                            if(bool(solver.Value(all_tasks_workers[product_id, station_id, workstep_id, worker_id].exists)) == True): 
                                assigned_products_workers[worker_id].append(
                                    assigned_task_type(
                                        start = solver.Value(all_tasks[product_id, station_id, workstep_id].start),
                                        exists = True,#exists = solver.Value(all_tasks[product_id, station_id, workstep_id].exists),
                                        product = product_id,
                                        station = station_id,
                                        workstep = workstep_id,
                                        duration = solver.Value(all_tasks[product_id, station_id, workstep_id].duration)))

            # Sort by starting time.
            for worker in all_workers:
                assigned_products_workers[worker].sort()
            for product in all_products:
                assigned_stations[product].sort()

            # Write the solution to the XML input file
            xml_writer = XML_Writer(config_path, init_path)
            xml_writer.write_solution(sol_path, assigned_products_workers, assigned_stations, solution_objective, "OR-Tools", solution_time, solution_status, 0)
            conf_f = config_path.split(".")[0].split("/")[4]
            ini_f = init_path.split(".")[0].split("/")[4]
            sol_f = sol_path.split(".")[0].split("/")[4]
            new_row = f"{conf_f},{ini_f},{sol_f},{solution_status},{solution_time},{solution_objective},{model_loading_time}, {horizon}\n"
            with open('Results.csv','a') as fd:
                fd.write(new_row)