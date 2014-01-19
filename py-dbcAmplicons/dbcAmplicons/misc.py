#### Misc functions

import os, errno

def make_sure_path_exists(path):
    if path != '':
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

