class XML_Writer(object):
    """Writes solution into XML Files"""

    def __init__(self, config_path, init_path):
        import xml.etree.ElementTree as etree
        self.config_path = config_path
        self.init_path = init_path
        self.configtree = etree.parse(self.config_path)
        self.configroot = self.configtree.getroot()
        self.inittree = etree.parse(self.init_path)
        self.initroot = self.inittree.getroot()

    def get_worker_dictionary(self):
        dict_workers = {}
        counter = 1
        for workers in self.initroot.findall("Workers"):
            for worker in workers.findall("Worker"):
                worker_id = worker.get("id")
                dict_update = {counter: worker_id}
                dict_workers.update(dict_update)
                counter += 1
        return(dict_workers)

    def get_variant_dictionary(self):
        dict_variants = {}
        counter = 1
        for variants in self.configroot.findall("Variants"):
            for variant in variants.findall("Variant"):
                variant_id = variant.get("id")
                dict_update = {counter: variant_id}
                dict_variants.update(dict_update)
                counter += 1
        return(dict_variants)

    def get_product_dicitionary(self):
        dict_products = {}
        counter = 1
        for products in self.initroot.findall("Products"):
            for product in products.findall("Product"):
                product_id = product.get("id")
                dict_update = {counter: product_id}
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
                    dict_update = {counter: station_id}
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
                    dict_update = {counter: buffer_id}
                    dict_buffers.update(dict_update)
                    counter += 1
        return(dict_buffers)

    def write_solution(self, sol_path, assigned_products_workers, assigned_stations, solution_objective, solver, solution_time, solution_status, solution_gap):
        import xml.etree.ElementTree as etree
        from datetime import datetime
        #from XML_Reader import XML_Reader
        def sortchildrenby(parent, attr):
            parent[:] = sorted(parent, key=lambda child: int(child.get(attr)))

        worker_dict = self.get_worker_dictionary()
        station_dict = self.get_station_dictionary()
        product_dict = self.get_product_dicitionary()
        variant_dict = self.get_variant_dictionary()
        buffer_dict = self.get_buffer_dictionary()


        tree = etree.ElementTree()
        root = etree.Element('FlowShopSolution', configuration='BoschExample', initialState='BoschExample-1')
        tree._setroot(root)

        events = etree.SubElement(root, 'Events')
        workers = etree.SubElement(events, 'Workers')
        products = etree.SubElement(events, 'Products')
        
        # Making the Worker Element in the Events Entry
        worker_events = []
        for i in range(len(assigned_products_workers)):
            counter = 0
            for j in assigned_products_workers[i]:
                if j.exists:
                    new_station = j.station
                    if counter == 0:
                        station = j.station
                        start_time = j.start
                        end_time = j.start + j.duration
                        counter += 1
                    elif station == new_station:
                        end_time = j.start + j.duration
                        counter += 1
                    else:
                        worker_events.append(etree.SubElement(workers, 'Event', worker=f'{worker_dict[i+1]}', station=f'{station_dict[station+1]}', start=f'{start_time}', end=f'{end_time}'))
                        start_time = j.start
                        end_time = j.start + j.duration
                        station = new_station
                        counter = 1
            worker_events.append(etree.SubElement(workers, 'Event', worker=f'{worker_dict[i+1]}', station=f'{station_dict[station+1]}', start=f'{start_time}', end=f'{end_time}'))

        sortchildrenby(workers, 'start')
        # Making the Product Element in the Events entry
        product_events = []
        for i in range(len(assigned_stations)):
            counter = 0
            for j in assigned_stations[i]:
                if j.exists:
                    new_station = j.station
                    if counter == 0:
                        station = j.station
                        start_time = j.start
                        end_time = j.start + j.duration
                        counter += 1
                    elif station == new_station:
                        end_time = j.start + j.duration
                        counter += 1
                    else:
                        product_events.append(etree.SubElement(products, 'Event', product=f'{product_dict[i+1]}', location=f'{station_dict[station+1]}', start=f'{start_time}', end=f'{end_time}'))
                        start_time = j.start
                        if start_time != end_time and station+1<len(station_dict):
                            product_events.append(etree.SubElement(products, 'Event', product=f'{product_dict[i+1]}', location=f'{buffer_dict[station+1]}', start=f'{end_time}', end=f'{start_time}'))
                        end_time = j.start + j.duration
                        station = new_station
                        counter = 1
            product_events.append(etree.SubElement(products, 'Event', product=f'{product_dict[i+1]}', location=f'{station_dict[station+1]}', start=f'{start_time}', end=f'{end_time}'))  
        sortchildrenby(products, 'start')

        # Making the Schedule Element
        schedule = etree.SubElement(root, 'Schedule', objective=f'{solution_objective}', solver=solver, status=f'{solution_status}', solvingTime=f'{solution_time}', gap=f'{solution_gap}')

        xml_reader = XML_Reader(self.config_path, self.init_path)
        workstep_dict = xml_reader.get_workstep_dictionary()
        variants = xml_reader.get_products_type()

        solution_worker = []
        solution_work_steps = []
        solution_work_step = []
        now = datetime.now()
        for worker in range(len(assigned_products_workers)):
            solution_worker.append(etree.SubElement(schedule, 'Worker', id=f'{worker_dict[worker+1]}'))
            solution_work_steps.append(etree.SubElement(solution_worker[worker], 'AssignedWorkSteps'))
            solution_work_step.append([])
            for w_step in assigned_products_workers[worker]:
                if w_step.exists:
                    stat = w_step.station + 1
                    prod = w_step.product + 1
                    workstep = w_step.workstep + 1
                    start_time = w_step.start
                    end_time = w_step.start + w_step.duration
                    workstep = [variants[prod-1]+1, stat]
                    workstep_id = 0
                    for wstep_id in workstep_dict:
                        if workstep_dict[wstep_id] == workstep:
                            workstep_id = wstep_id
                    solution_work_step[worker].append(etree.SubElement(solution_work_steps[worker], 'Workstep', workstep=f'{workstep_id}', station=f'{station_dict[stat]}', product=f'{product_dict[prod]}', 
                                                                       start=f'{start_time}', end=f'{end_time}'))
        
        tree.write(sol_path)