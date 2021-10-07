"""Logging functions for FMI-PPN"""

import logging

_logger = None

_logger_severity = {
    'critical': logging.CRITICAL,
    'debug': logging.DEBUG,
    'error': logging.ERROR,
    'info': logging.INFO,
    'warning': logging.WARNING,
}

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
    if _logger is None:
        raise RuntimeError("Tried to write to log before logging was configured! "
                           "Call 'configure_logging' first.")
    lvl = level.lower()
    severity = _logger_severity.get(lvl)
    if severity is None:
        raise ValueError("Unknown logging level {}".format(level))
    _logger.log(severity, msg, *args, **kwargs)
