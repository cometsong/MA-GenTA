import string
import datetime

import logbook
from logbook import Logger, FileHandler
from logbook.more import ColorizedStderrHandler

logbook.set_datetime_format('local')

log_name = 'Targeted_Pipeline'
show_level = 'INFO'
log_level = 'DEBUG'

# DEFAULT_FORMAT_STRING = ('[{record.time:%Y-%m-%d %H:%M:%S.%f%z}] '
#                          '{record.level_name}: {record.channel}: {record.message}')
FORMAT_STRING = ('[{record.time:%Y%m%d-%H:%M:%S}]'
                 ' {record.level_name:^7}:' # centered within 7
                 ' {record.channel}:' # channel = log.name
                 ' {record.message}'
                 )

def inject_extra_group(record):
    """Add name of "subgroup" into 'extras' dict"""
    record.extra['group'] = ''


def log_file_init(log_name=__name__, logfile=None):
    NOW = datetime.datetime.now().isoformat('T', 'seconds').replace(':','')
    if not logfile:
        _puncts = string.punctuation.replace('_-','') # allow specific punctuation in Title
        logbase = log_name.replace(_puncts,'').replace(' ','')
        logfile = '.'.join([NOW, logbase, 'log'])
    return logfile


def log_init(name=__name__, level='NOTICE', show_level=None,
             format_string=FORMAT_STRING, logfile=None):
    """Initialize a new Logger to file and colorized stderr stream"""
    logfile = log_file_init(log_name=name, logfile=logfile)

    file_handler = FileHandler(logfile, level=level, format_string=format_string, bubble=True)
    show_level = show_level if show_level else level
    cstd_handler = ColorizedStderrHandler(level=show_level, format_string=format_string, bubble=False)

    level = logbook.lookup_level(level)

    logger = Logger(name, level=level)
    logger.handlers.append(file_handler)
    logger.handlers.append(cstd_handler)

    return logger


log = log_init(name=log_name, level=log_level, show_level=show_level)

