class ORSolver(object):
    """Solver Class of the JobShop Problem"""

    def solve_with_disunctive_model(self, config_path, init_path, sol_path, max_time):
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
      print('Loading data files')
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
      print(f'  #Stations: {stations_count}')
      all_stations = range(stations_count)

      products_count = len(products_type)
      print(f'  #Products: {products_count}')
      all_products = range(products_count)

      num_product_types = 4  # len(production_time_matrix)
      print(f'  #Product types: {num_product_types}')

      horizon = 1000000
      print(f'  Horizon: {horizon}')

      workers_count = len(worker_skill)
      print(f'  #Workers: {workers_count}')
      all_workers = range(workers_count)

      # Named tuple to store information about created variables.
      worker_task_type = collections.namedtuple(
          'worker_task_type', 'exists start end duration station product step')

      task_type = collections.namedtuple(
          'task_type', 'start end exists_vars worker_durations')

      # Named tuple to manipulate solution information.
      assigned_task_type = collections.namedtuple(
          'assigned_task_type', 'start exists product station workstep duration')

      # Creates job intervals, for two intervals assigned to the same worker,
      # whether one occurs before the other(is_before), whether both intervals
      # exist for the same worker (both_exist) and the start_time_diff (the
      # difference in start times between two jobs if they both exist for this
      # worker, otherwise we set this to 0) and add to the corresponding machine
      # lists.
      print('\nBuilding model')
      all_tasks = {}  # indexed by (product, station, workstep)
      all_worker_tasks = collections.defaultdict(list)  # indexed by workers
      task_starts = []
      print(f'  Task variables: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station in range(products_initial[0][product], stations_count):
          for workstep in range(
              number_manual_worksteps[products_type[product]][station]):
            # Initial state: product is performed until products_initial[1][product]
            # and the first step to do is products_initial[0][product]
            if station == products_initial[0][product] and (
                workstep < products_initial[1][product] - 1):
              continue
            suffix = '_%i_%i_%i' % (product, station, workstep)
            # start time of each job
            start_var = model.NewIntVar(0, horizon, 'start' + suffix)
            # end time of each job
            end_var = model.NewIntVar(0, horizon, 'end' + suffix)

            # dictionary of all tasks
            current_task = task_type(
                start=start_var,
                end=end_var,
                exists_vars=[],
                worker_durations=[])
            all_tasks[product, station, workstep] = current_task
            task_starts.append(start_var)
            for worker in all_workers:
              # for each worker we create an optional interval for each job
              suffix_worker = '_%i_%i_%i_%i' % (worker, product, station, workstep)
              # variable saying whether the given job is produced by the worker
              exists_var = model.NewBoolVar('exists' + suffix_worker)
              # duration for the jobs
              duration_worker = int(
                  round(
                      float(production_time_matrix[products_type[product]][station]
                            [workstep]) * float(100) / float(worker_skill[worker])))
              # dictionary of optional intervals
              all_worker_tasks[worker].append(
                  worker_task_type(
                      exists=exists_var,
                      start=start_var,
                      end=end_var,
                      duration=duration_worker,
                      station=station,
                      product=product,
                      step=workstep))
              # Append worker info on the current task.
              current_task.exists_vars.append(exists_var)
              current_task.worker_durations.append(duration_worker)

      tasks = range(len(all_worker_tasks[0]))

      # end = start + duration
      print(f'  Task durations: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station in range(products_initial[0][product], stations_count):
          for workstep in range(
              number_manual_worksteps[products_type[product]][station]):
            if station == products_initial[0][product] and (
                workstep < products_initial[1][product] - 1):
              continue
            task = all_tasks[product, station, workstep]
            min_duration = min(task.worker_durations)
            shifted = [d - min_duration for d in task.worker_durations]
            model.Add(task.end == task.start + min_duration +
                      cp_model.LinearExpr.ScalProd(task.exists_vars, shifted))

      # each interval is present for exactly one worker
      print(f'  One worker per task: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station in range(products_initial[0][product], stations_count):
          for workstep in range(
              number_manual_worksteps[products_type[product]][station]):
            if station == products_initial[0][product] and (
                workstep < products_initial[1][product] - 1):
              continue
            model.Add(sum(all_tasks[product, station, workstep].exists_vars) == 1)

      # if the release time for a worker is greater than 0, then he can't start
      # working until then
      print(f'  Release time: {time.perf_counter() - tic:0.4f}s')
      for worker in all_workers:
        if releasetime_workers[worker] > 0:
          for task in tasks:
            model.Add(all_worker_tasks[worker][task].start >=
                      releasetime_workers[worker] + worker_movement_time_matrix[
                          workers_initial_station[worker]][station]).OnlyEnforceIf(
                              all_worker_tasks[worker][task].exists)

      # if a product has some leftover time at a station, the product cannot start
      # being processed elsewhere before the releasetime and no other product can be
      # produced at that station until the product is done
      print(f'  Leftover constraints: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station in range(products_initial[0][product], stations_count):
          if leftover_job_products[0][product] <= 0:
            continue
          for workstep in range(
              number_manual_worksteps[products_type[product]][station]):
            if station == products_initial[0][product] and (
                workstep < products_initial[1][product] - 1):
              continue
            model.Add(
                all_tasks[product, station,
                          workstep].start >= leftover_job_products[0][product])
            # the station is blocked for all further products until the release time
            for product2 in range(product + 1, products_count):
              model.Add(all_tasks[product2, leftover_job_products[1][product] - 1,
                                  0].start >= leftover_job_products[0][product])
            # the next buffer must be clear when the product finishes (but only if
            # the final workstep at the station is being completed before the
            # release time)
            if products_initial[1] == 1:
              for product1 in range(0, product):
                if products_initial[0][product1] == leftover_job_products[1][
                    product]:
                  continue
                model.Add(all_tasks[product1, leftover_job_products[1][product],
                                    0].start <= leftover_job_products[0][product])

      # precedence constraints stations
      print(f'  Station precedendce constraints: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station1 in range(products_initial[0][product], stations_count - 1):
          station2 = station1 + 1  # for station2 in range(station1+1,stations_count):
          model.Add(
              all_tasks[product, station1,
                        number_manual_worksteps[products_type[product]][station1] -
                        1].end <= all_tasks[product, station2, 0].start)

      # precedence constraints products
      print(f'  Product precedence constraints {time.perf_counter() - tic:0.4f}s')
      for product1 in range(0, products_count - 1):  # all_products:
        product2 = product1 + 1  # for product2 in range(product1+1, products_count):
        for station in range(products_initial[0][product1], stations_count):
          if (station == products_initial[0][product2] and
              products_initial[1][product2] != 1):
            continue
          model.Add(
              all_tasks[product1, station,
                        number_manual_worksteps[products_type[product1]][station] -
                        1].end <= all_tasks[product2, station, 0].start)

      # buffer constraint
      print(f'  Buffer constraints: {time.perf_counter() - tic:0.4f}s')
      for product1 in range(0, products_count - 1):  # all_products:
        product2 = product1 + 1  # for product2 in range(product1+1, products_count):
        for station in range(
            max([products_initial[0][product1] - 1, products_initial[0][product2]]),
            stations_count - 1):
          if (station + 1 == products_initial[0][product1] and
              products_initial[1][product1] != 1):
            continue
          model.Add(all_tasks[product1, station + 1, 0].start <= all_tasks[
              product2, station,
              number_manual_worksteps[products_type[product2]][station] - 1].end)

      # automatic station
      print(f'  Automatic stations: {time.perf_counter() - tic:0.4f}s')
      for product in all_products:
        for station in range(products_initial[0][product], stations_count):
          if (station == products_initial[0][product] and
              0 < products_initial[1][product] - 1):
            continue
          if number_manual_worksteps[products_type[product]][station] == 2:
            model.Add(
                all_tasks[product, station,
                          0].end == all_tasks[product, station, 1].start -
                automatic_worksteps_times_matrix[products_type[product]][station])

      # Disjunctions between tasks.
      print(f'  Time disjunctions: {time.perf_counter() - tic:0.4f}s')
      constraint_counter = 0
      for worker in all_workers:
        for task1 in tasks:
          for task2 in range(task1 + 1, len(all_worker_tasks[worker])):
            suffix = '_%i_%i_%i' % (worker, task1, task2)
            station1 = all_worker_tasks[worker][task1].station
            station2 = all_worker_tasks[worker][task2].station
            product1 = all_worker_tasks[worker][task1].product
            product2 = all_worker_tasks[worker][task2].product
            if ((station1 >= station2) and (product1 >= product2)) and not (
                (station1 == station2 + 1) and
                (product1 == product2)) and not ((station1 == station2) and
                                                 (product1 == product2 + 1)):
              continue
            if ((station1 <= station2) and (product1 <= product2)) and not (
                (station1 == station2 - 1) and
                (product1 == product2)) and not ((station1 == station2) and
                                                 (product1 == product2 - 1)):
              continue
            constraint_counter += 1
            is_before_var = model.NewBoolVar('is_before' + suffix)
            exists1 = all_worker_tasks[worker][task1].exists
            exists2 = all_worker_tasks[worker][task2].exists
            model.Add(
                all_worker_tasks[worker][task2].start >=
                all_worker_tasks[worker][task1].end +
                worker_movement_time_matrix[station1][station2]).OnlyEnforceIf(
                    [is_before_var, exists1, exists2])
            model.Add(all_worker_tasks[worker][task1].start >= (
                all_worker_tasks[worker][task2].end +
                worker_movement_time_matrix[station2][station1])).OnlyEnforceIf(
                    [is_before_var.Not(), exists1, exists2])
      print('    #is_before constraints: ', constraint_counter)

      # Initial positions of the workers means that no production step can start
      # before the worker has walked to the station from its initial position.
      print(f'  Initial worker positions: {time.perf_counter() - tic:0.4f}s')
      for worker in all_workers:
        for task in tasks:
          model.Add(all_worker_tasks[worker][task].start >=
                    worker_movement_time_matrix[workers_initial_station[worker]][
                        all_worker_tasks[worker][task].station]).OnlyEnforceIf(
                            all_worker_tasks[worker][task].exists)

      # Makespan objective.
      print(f'  Objective: {time.perf_counter() - tic:0.4f}s')
      obj_var = model.NewIntVar(0, horizon, 'makespan')
      model.Add(obj_var == all_tasks[
          products_count - 1, stations_count - 1,
          number_manual_worksteps[products_type[products_count -
                                                1]][stations_count - 1] - 1].end)
      model.Minimize(obj_var)

      # Solve model.
      tada = time.perf_counter()
      model_loading_time = tada - tic
      print('\nSolve model')
      solver = cp_model.CpSolver()
      solver.parameters.num_search_workers = 16  # if multiple workers
      solver.parameters.max_time_in_seconds = max_time - model_loading_time
      solver.parameters.log_search_progress = True
      status = solver.Solve(model)
      toc = time.perf_counter()

      if status == cp_model.INFEASIBLE:
        print(f'Shown infeasible in {toc - tic:0.4f} seconds')
        status_message = 'INFEASIBLE'
        solution_objective = horizon
        solution_time = toc - tic
        solution_status = 'INFEASIBLE'

      if status == cp_model.MODEL_INVALID:
        print(f'Shown model invalid in {toc - tic:0.4f} seconds')
        status_message = 'MODEL_INVALID'
        status_message = 'UNSOLVED'
        solution_objective = horizon
        solution_time = toc - tic
        solution_status = 'UNSOLVED'

      if status == cp_model.OPTIMAL:
        print(f'Solved optimally in {toc - tic:0.4f} seconds')
        status_message = 'OPTIMAL'
        # Create one list of assigned tasks per machine.
      if status == cp_model.FEASIBLE:
        print(f'Solved suboptimally in {toc - tic:0.4f} seconds')
        status_message = 'SUBOPTIMAL'

      if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        # Finally print the solution found.
        solution_objective = solver.ObjectiveValue()
        solution_time = toc - tic
        solution_status = status_message
        #print('Schedule Length: %i' % solution_objective)

        assigned_stations = collections.defaultdict(list)
        for product_id in all_products:
          for station_id in all_stations:
            if station_id < products_initial[0][product_id]:
              continue
            for workstep_id in range(
                number_manual_worksteps[products_type[product_id]][station_id]):
              if station_id == products_initial[0][product_id] and (
                  workstep_id < products_initial[1][product_id] - 1):
                continue
              task = all_tasks[product_id, station_id, workstep_id]

              d = -1
              for var, duration in zip(task.exists_vars, task.worker_durations):
                if solver.BooleanValue(var):
                  d = duration

              assigned_stations[product_id].append(
                  assigned_task_type(
                      start=solver.Value(task.start),
                      exists=True,  # exists=solver.Value(all_tasks[product_id, station_id, workstep_id].exists),
                      product=product_id,
                      station=station_id,
                      workstep=workstep_id,
                      duration=d))

        # Create one list of assigned tasks per worker
        assigned_products_workers = collections.defaultdict(list)
        for worker_id in all_workers:
          for task in tasks:
            if solver.BooleanValue(all_worker_tasks[worker_id][task].exists):
              assigned_products_workers[worker_id].append(
                  assigned_task_type(
                      start=solver.Value(all_worker_tasks[worker_id][task].start),
                      exists=True,
                      product=all_worker_tasks[worker_id][task].product,
                      station=all_worker_tasks[worker_id][task].station,
                      workstep=all_worker_tasks[worker_id][task].step,
                      duration=solver.Value(
                          all_worker_tasks[worker_id][task].duration)))
        # Sort by starting time.
        for worker in all_workers:
          assigned_products_workers[worker].sort()
        for product in all_products:
          assigned_stations[product].sort()

        # Write the solution to the XML input file
        xml_writer = XML_Writer(config_path, init_path)
        xml_writer.write_solution(sol_path, assigned_products_workers, assigned_stations, solution_objective, "OR-Tools", solution_time, solution_status, 0)
      else:
        solution_objective = horizon
        solution_time = toc - tic
        solution_status = 'UNSOLVED'
      conf_f = config_path.split(".")[0].split("/")[4]
      ini_f = init_path.split(".")[0].split("/")[4]
      sol_f = sol_path.split(".")[0].split("/")[4]
      new_row = f"{conf_f},{ini_f},{sol_f},{solution_status},{solution_time},{solution_objective},{model_loading_time}, {horizon}\n"
      with open('Results.csv','a') as fd:
            fd.write(new_row)
