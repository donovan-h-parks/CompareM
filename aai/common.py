###############################################################################
#
# common.py - utility functions used in many places in CheckM
#
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

import os
import errno
import sys
import logging


def check_file_exists(input_file):
    """Check if file exists."""
    if not os.path.exists(input_file) and os.path.isfile(input_file):
        logger = logging.getLogger()
        logger.error('  [Error] Input file does not exists: ' + input_file + '\n')
        sys.exit()


def check_dir_exists(input_dir):
    """Check if directory exists."""
    if not os.path.exists(input_dir) and os.path.isdir(input_dir):
        logger = logging.getLogger()
        logger.error('  [Error] Input directory does not exists: ' + input_dir + '\n')
        sys.exit()


def make_sure_path_exists(path):
    """Create directory if it does not exist."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger()
            logger.error('  [Error] Specified path does not exist: ' + path + '\n')
            sys.exit()
