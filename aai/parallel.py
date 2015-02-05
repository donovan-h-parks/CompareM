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

    def __producer(self, producerCallback, producerQueue, consumerQueue):
        """Process data items with worker thread."""
        while True:
            dataItem = producerQueue.get(block=True, timeout=None)
            if dataItem == None:
                break

            rtn = producerCallback(dataItem)

            consumerQueue.put(rtn)

    def __processManager(self, cpus, producerCallback, producerQueue, consumerQueue):
        """Setup producer processes and manage completion of processes."""

        try:
            producerProc = [mp.Process(target=self.__producer, args=(producerCallback, producerQueue, consumerQueue)) for _ in range(cpus)]

            for p in producerProc:
                p.start()

            for p in producerProc:
                p.join()

            consumerQueue.put(None)
        except:
            for p in producerProc:
                p.terminate()

    def run(self, producer, consumer, dataItems, cpus):
        """Setup producers and consumer processes."""

        # populate producer queue with data to process
        consumerQueue = mp.Queue()
        producerQueue = mp.Queue()
        for d in dataItems:
            producerQueue.put(d)

        for _ in range(cpus):
            producerQueue.put(None)  # signal processes to terminate

        try:

            mangerProc = mp.Process(target=self.__processManager, args=(cpus, producer, producerQueue, consumerQueue))

            mangerProc.start()

            # process items produced by workers
            processedItems = 0
            while True:
                rtn = consumerQueue.get(block=True, timeout=None)
                if rtn == None:
                    break

                consumer(rtn)

                processedItems += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) items.' % (processedItems, len(dataItems), float(processedItems) * 100 / len(dataItems))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            sys.stdout.write('\n')

            mangerProc.join()
        except:
            mangerProc.terminate()
