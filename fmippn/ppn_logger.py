"""Logging functions for FMI-PPN"""

import logging

_logger = None

def config_logging(fname, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S"):
    """Setup logger object using logging.basicConfig.

    Input:
        fname -- name of the log file (str)

    `level` and `datefmt` arguments correspond to logging.basicConfig arguments.
    """
    global _logger
    formatter = logging.Formatter(
        fmt="%(asctime)s (%(name)s) %(levelname)s: %(message)s",
        datefmt=datefmt
    )
    handler = logging.FileHandler(fname)
    handler.setFormatter(formatter)
    _logger = logging.getLogger("FMIPPN")
    _logger.addHandler(handler)
    _logger.setLevel(level)

def write_to_log(level, msg, *args, **kwargs):
    """Write `msg` at logging level `level` to log.
    args and kwargs are passed to logging functions.
    """
    lvl = level.lower()
    if lvl == 'critical':
        _logger.critical(msg, *args, **kwargs)
    elif lvl == 'error':
        _logger.error(msg, *args, **kwargs)
    elif lvl == 'warning':
        _logger.warning(msg, *args, **kwargs)
    elif lvl == 'info':
        _logger.info(msg, *args, **kwargs)
    elif lvl == 'debug':
        _logger.debug(msg, *args, **kwargs)
    else:
        raise ValueError("Unknown logging level {}".format(level))
