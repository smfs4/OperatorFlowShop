class XML_Reader(object):
    """Provides Methods for extracting necessary information from the XML file"""
    def __init__(self, config_path, init_path):
        import xml.etree.ElementTree as ET
        self.configfilepath = config_path
        self.initfilepath = init_path
        self.configtree = ET.parse(self.configfilepath)
        self.configroot = self.configtree.getroot()
        self.inittree = ET.parse(self.initfilepath)
        self.initroot = self.inittree.getroot()

    def get_num_workers(self):
        num_workers = len(self.initroot.findall("Workers")[0])
        return(num_workers)

    def get_num_products(self):
        num_products = len(self.initroot.findall("Products")[0])
        return(num_products)

    def get_num_variants(self):
        num_variants = len(self.configroot.findall("Variants")[0])
        return(num_variants)

    def get_num_stations(self):
        num_stations = 0
        for i in self.configroot.findall("ProcessingUnits")[0]:
            if i.get('type') == "STATION":
                num_stations += 1
        return(num_stations)

    def get_worker_dictionary(self):
        dict_workers = {}
        counter = 1
        for workers in self.initroot.findall("Workers"):
            for worker in workers.findall("Worker"):
                worker_id = worker.get("id")
                if worker_id not in dict_workers:
                    dict_update = {worker_id: counter}
                    dict_workers.update(dict_update)
                    counter += 1
        return(dict_workers)

    def get_variant_dictionary(self):
        dict_variants = {}
        counter = 1
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                variant_id = variant.get("id")
                if variant_id not in dict_variants:
                    dict_update = {variant_id: counter}
                    dict_variants.update(dict_update)
                    counter += 1
        return(dict_variants)

    def get_product_dicitionary(self):
        dict_products = {}
        counter = 1
        for products in self.initroot.findall("Products"):
            for product in products.findall("Product"):
                product_id = product.get("id")
                if product_id not in dict_products:
                    dict_update = {product_id: counter}
                    dict_products.update(dict_update)
                    counter += 1
        return(dict_products)

    def get_station_dictionary(self):
        dict_stations = {}
        counter = 1
        for processing_units in self.configroot.findall("ProcessingUnits"):
            for processing_unit in processing_units.findall("ProcessingUnit"):
                if processing_unit.get("type") == "STATION":
                    station_id = processing_unit.get("id")
                    if station_id not in dict_stations:
                        dict_update = {station_id: counter}
                        dict_stations.update(dict_update)
                        counter += 1
        return(dict_stations)

    def get_buffer_dictionary(self):
        dict_buffers = {}
        counter = 1
        for processing_units in self.configroot.findall("ProcessingUnits"):
            for processing_unit in processing_units.findall("ProcessingUnit"):
                if processing_unit.get("type") == "BUFFER":
                    buffer_id = processing_unit.get("id")
                    if buffer_id not in dict_buffers:
                        dict_update = {buffer_id: counter}
                        dict_buffers.update(dict_update)
                        counter += 1
        return(dict_buffers)

    def get_workstep_dictionary(self):
        dict_workstep = {}
        dict_stations = self.get_station_dictionary()
        dict_variant = self.get_variant_dictionary()
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                for worksteps in variant.findall("WorkSteps"):
                    for workstep in worksteps.findall("WorkStep"):
                        location = workstep.get("location")
                        workstep_id = workstep.get("id")
                        automatic = workstep.get("automatic")
                        variant_no = dict_variant[variant.get("id")]
                        station_no = dict_stations[location]
                        dict_update = {workstep_id:[variant_no, station_no]}
                        dict_workstep.update(dict_update)
        return(dict_workstep)



    def get_worker_skill(self):
        worker_skill = []
        for workers in self.initroot.findall("Workers"):
            for worker in workers.findall("Worker"):
                worker_skill.append(int(worker.get("skill")))
        return(worker_skill)


    def get_worker_initial_stat(self):
        workers_initial_station = []
        station_dict = self.get_station_dictionary()
        for workers in self.initroot.findall("Workers"):
            for worker in workers.findall("Worker"):
                in_station = worker.get("location")
                workers_initial_station.append(station_dict[in_station]-1)
        return(workers_initial_station)

    def get_products_type(self):
        products_type = []
        variant_dict = self.get_variant_dictionary()
        for products in self.initroot.findall("Products"):
            for product in products.findall("Product"):
                var = product.attrib.get("variant")
                products_type.append(variant_dict[var]-1)
        return(products_type)

    def get_products_initial(self):
        import numpy as np
        num_products = self.get_num_products()
        products_initial = np.zeros((2, num_products)).astype(int)
        dict_stations = self.get_station_dictionary()
        dict_buffers = self.get_buffer_dictionary()
        dict_products = self.get_product_dicitionary()

        # first workstep to be scheduled as station and workstep for each product
        for products in self.initroot.findall("Products"):
            for product in products.findall("Product"):
                prod = product.get("id")
                prod_id = dict_products[prod]
                if product.get("releasetime") != None:
                    if product.get("step") == "manual" or product.get("step") == "takedown":
                        products_initial[0, prod_id-1] = dict_stations[product.get("location")]
                        products_initial[1, prod_id-1] = 1
                    else:
                        products_initial[0, prod_id-1] = dict_stations[product.get("location")]-1
                        products_initial[1, prod_id-1] = 2
                elif product.get("location") != None:
                    products_initial[0, prod_id-1] = dict_buffers[product.get("location")]
                    products_initial[1, prod_id-1] = 1
                else:
                    products_initial[0, prod_id-1] = 0
                    products_initial[1, prod_id-1] = 1
        products_initial = products_initial.tolist()
        return(products_initial)

    def get_production_time_matrix(self):
        import numpy as np
        maxWorksteps = self.get_max_worksteps()
        numVariants = self.get_num_variants()
        numStations = self.get_num_stations()
        dict_stations = self.get_station_dictionary()
        production_time_matrix = np.zeros((numVariants, numStations, maxWorksteps)).astype(int)
        var = 0
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                wstep = 0
                for worksteps in variant.findall("WorkSteps"):
                    for workstep in worksteps.findall("WorkStep"):
                        wstep += 1
                        automatic = workstep.attrib.get('automatic')
                        station = workstep.get("location")
                        stat = dict_stations[station]
                        for duration in workstep.findall("Duration"):
                            if automatic == "True":
                                for setuptime in workstep.findall("SetupTime"):
                                    production_time_matrix[var][stat-1][0] = int(setuptime.text)
                                for takedown_time in workstep.findall("TakedownTime"):
                                    production_time_matrix[var][stat-1][1] = int(takedown_time.text)
                            else:
                                production_time_matrix[var][stat-1][0] = int(duration.text)
                var += 1
        production_time_matrix = production_time_matrix.tolist()
        return(production_time_matrix)

    def get_automatic_times_matrix(self):
        import numpy as np
        numVariants = self.get_num_variants()
        numStations = self.get_num_stations()
        automatic_worksteps_times_matrix = np.zeros(shape=(numVariants,numStations)).astype(int)
        dict_stations = self.get_station_dictionary()
        var = 0
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                wstep = 0
                for worksteps in variant.findall("WorkSteps"):
                    for workstep in worksteps.findall("WorkStep"):
                        wstep += 1
                        automatic = workstep.attrib.get('automatic')
                        station = workstep.get("location")
                        if automatic == "True":
                            for duration in workstep.findall("Duration"):
                                stat = dict_stations[station]
                                automatic_worksteps_times_matrix[var][stat-1] = int(duration.text)
                var += 1
        automatic_worksteps_times_matrix = automatic_worksteps_times_matrix.tolist()
        return(automatic_worksteps_times_matrix)

    def get_max_worksteps(self):
        number_worksteps = self.get_number_worksteps()
        maxWorksteps = max(max(number_worksteps))
        return(maxWorksteps)

    def get_number_worksteps(self):
        import numpy as np
        numVariants = self.get_num_variants()
        numStations = self.get_num_stations()
        number_manual_worksteps = np.zeros(shape=(numVariants,numStations)).astype(int)
        dict_stations = self.get_station_dictionary()
        var = 0
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                wstep = 0
                for worksteps in variant.findall("WorkSteps"):
                    for workstep in worksteps.findall("WorkStep"):
                        wstep += 1
                        automatic = workstep.attrib.get('automatic')
                        station = workstep.get("location")
                        stat = dict_stations[station]
                        if automatic == "True":
                            number_manual_worksteps[var][stat-1] = 2
                        else:
                            number_manual_worksteps[var][stat-1] = 1
                var += 1
        number_manual_worksteps = number_manual_worksteps.tolist()
        return(number_manual_worksteps)

    def get_releasetime_workers(self):
        import numpy as np
        numWorkers = self.get_num_workers()
        dict_workers = self.get_worker_dictionary()
        releasetime = np.zeros((numWorkers)).astype(int)
        for workers in self.initroot.findall('Workers'):
            for worker in workers.findall('Worker'):
                work = worker.get("id")
                worker_id = dict_workers[work]
                if worker.get("releasetime") != None:
                    rem_time = int(worker.get("releasetime"))
                    releasetime[worker_id-1] = rem_time
        releasetime = releasetime.tolist()
        return(releasetime)

    def get_leftover_jobs_products(self):
        import numpy as np
        numProducts = self.get_num_products()
        dict_products = self.get_product_dicitionary()
        dict_stations = self.get_station_dictionary()
        leftover_job_products = np.zeros((2, numProducts)).astype(int)
        for products in self.initroot.findall("Products"):
            for product in products.findall("Product"):
                station = 0
                rel_time = 0
                if product.get("releasetime") != None:
                    rel_time = int(product.get("releasetime"))
                    station = int(dict_stations[product.get("location")])
                prod = product.get("id")
                product_id = dict_products[prod]
                leftover_job_products[0, product_id-1] = rel_time
                leftover_job_products[1, product_id-1] = station
        leftover_job_products = leftover_job_products.tolist()
        return(leftover_job_products)

    def get_worker_movement_matrix(self):
        import numpy as np
        numStations = self.get_num_stations()
        dict_stations = self.get_station_dictionary()
        worker_movement_time_matrix = np.zeros(shape=(numStations,numStations)).astype(int)
        for workerTT in self.configroot.findall('WalkingTimes'):
            for child in workerTT:
                travelTime = int(child.get("traveltime"))
                startS = child.get('start')
                startS = dict_stations[startS]
                endS = child.get('end')
                endS = dict_stations[endS]
                worker_movement_time_matrix[startS-1, endS-1] = travelTime

        worker_movement_time_matrix = worker_movement_time_matrix.tolist()
        return(worker_movement_time_matrix)