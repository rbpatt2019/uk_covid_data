# -*- coding: utf-8 -*-
"""Configure logger."""
import logging
import sys


class RedirectStream:
    """Redirect stdout and stderr to a logger file.

    From https://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    """

    def __init__(self, logger: logging.Logger, level: int) -> None:
        """Initialise the RedirectStream object.

        Parameters
        ----------
        logger : logging.Logger
            A logger with a FileHandler
        level : int
            The logging level to associate with the re-direction
        """
        self.logger = logger
        self.level = level
        self.linebuf = ""

    def write(self, buf: str) -> None:
        """Write a redirection to the file.

        Parameters
        ----------
        buf : str
            The message to write to the file
        """
        self.logger.log(self.level, buf.rstrip())

    def flush(self) -> None:
        """Meet structure requirements.

        Serves no purpose other than to meet some expectations of certain packages.
        """
        pass


def get_logger(module: str, file: str, redirect: bool) -> logging.Logger:
    """Configure a file logger for use in a script.

    Parameters
    ----------
    module : str
        The name of the module from which the logger is called
    file : str
        The name of the log file to which the logger will write
    redirect : bool
        Whether to redirect sys.stdout and sys.stderr to the file.

    Returns
    -------
    logging.Logger
        The configured logger instance.
    """
    logger = logging.getLogger(module)

    handler = logging.FileHandler(file)
    formatter = logging.Formatter(
        "{asctime} :: {levelname} :: {name} :: {message}", style="{"
    )

    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)

    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    if redirect:
        sys.stdout = RedirectStream(logger, logging.INFO)  # type: ignore[assignment]
        sys.stderr = RedirectStream(logger, logging.ERROR)  # type: ignore[assignment]

    return logger
