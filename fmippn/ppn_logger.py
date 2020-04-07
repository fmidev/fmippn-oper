"""Logging functions for FMI-PPN"""

import logging

def config_logging(fname, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"):
    """Setup logger object using logging.basicConfig.

    Input:
        fname -- name of the log file (str)

    `level` and `datefmt` arguments correspond to logging.basicConfig arguments.
    """
    logging.basicConfig(
        filename=fname,
        level=level,
        datefmt=datefmt,
        format="%(asctime)s (%(name)s) %(levelname)s: %(message)s",
    )

def write_to_log(level, msg, *args, **kwargs):
    """Write `msg` at logging level `level` to log.
    args and kwargs are passed to logging functions.
    """
    lvl = level.lower()
    if lvl == 'critical':
        logging.critical(msg, *args, **kwargs)
    elif lvl == 'error':
        logging.error(msg, *args, **kwargs)
    elif lvl == 'warning':
        logging.warning(msg, *args, **kwargs)
    elif lvl == 'info':
        logging.info(msg, *args, **kwargs)
    elif lvl == 'debug':
        logging.debug(msg, *args, **kwargs)
    else:
        raise ValueError("Unknown logging level {}".format(level))
