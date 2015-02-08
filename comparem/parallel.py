###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import sys
import logging
import multiprocessing as mp


class Parallel(object):
    """Processes data in parallel.

    This class implements a producer/consumer model for
    processing data in parallel. In order to keep data
    handling efficient the main process acts as the
    consumer. This keeps all consumed data in the
    address space of the main process.

    Producer and consumer functions are passed to
    the run() method. The return value of the producer
    function is provided to the consumer as an argument.

    Example:

        class Test(object):
            def __init__(self):
                self.total = 0

            def _producer(self, data):
                return data * data

            def _consumer(self, dataSquard):
                self.total += dataSquard

            def run(self):
                parallel = Parallel()
                parallel.run(self._producer, self._consumer, [1, 2, 3, 4, 5], 2)

                print self.total
    """

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def __producer(self, producer_callback, producer_queue, consumer_queue):
        """Process data items with producer processes.

        Parameters
        ----------
        producer_callback : function
            Function to process data items.
        producer_queue : queue
            Data items to process.
        consumer_queue : queue
            Queue for holding processed data items to be consumed.
        """
        while True:
            dataItem = producer_queue.get(block=True, timeout=None)
            if dataItem == None:
                break

            rtn = producer_callback(dataItem)

            consumer_queue.put(rtn)

    def __process_manager(self, cpus, producer_callback, producer_queue, consumer_queue):
        """Setup producer processes and manage completion of processes.

        Parameters
        ----------
        cpus : int
            Number of processes to create.
        producer_callback : function
            Function to process data items.
        producer_queue : queue
            Data items to process.
        consumer_queue : queue
            Queue for holding processed data items to be consumed.

        """

        try:
            producer_proc = [mp.Process(target=self.__producer, args=(producer_callback, producer_queue, consumer_queue)) for _ in range(cpus)]

            for p in producer_proc:
                p.start()

            for p in producer_proc:
                p.join()

            consumer_queue.put(None)
        except:
            for p in producer_proc:
                p.terminate()

    def run(self, producer, consumer, data_items, cpus):
        """Setup producers and consumer processes.


        Parameters
        ----------
        producer : function
            Function to process data items.
        consumer : queue
            Function to consumed processed data items.
        data_items : list
            Data items to process.
        cpus : int
            Number of processes to create.
        """

        # populate producer queue with data to process
        producer_queue = mp.Queue()
        for d in data_items:
            producer_queue.put(d)

        for _ in range(cpus):
            producer_queue.put(None)  # signal processes to terminate

        try:
            consumer_queue = mp.Queue()
            manager_proc = mp.Process(target=self.__process_manager, args=(cpus, producer, producer_queue, consumer_queue))

            manager_proc.start()

            # process items produced by workers
            processed_items = 0
            while True:
                status = '    Finished processing %d of %d (%.2f%%) items.' % (processed_items, len(data_items), float(processed_items) * 100 / len(data_items))
                sys.stdout.write('%s\r' % status)
                sys.stdout.flush()

                rtn = consumer_queue.get(block=True, timeout=None)
                if rtn == None:
                    break

                consumer(rtn)

                processed_items += 1

            sys.stdout.write('\n')

            manager_proc.join()
        except:
            self.logger.warning('  [Warning] Exception encountered while processing data.')
            manager_proc.terminate()
