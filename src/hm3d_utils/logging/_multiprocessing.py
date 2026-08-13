
import multiprocessing
import logging 
from typing import Tuple
from logging.handlers import QueueHandler


def create_shared_logger_process_for_multiprocessing(queue: multiprocessing.Queue, log_file_path: str) -> None:
    """
    Function that runs in the logger process, continuously retrieving log records from the 
    shared logging queue and handling them by writing to both a log file and the console. 
    This function processes log records passed through the queue from worker processes, 
    ensuring logs are stored in a file and displayed on the console in real time.

    Args:
        queue (multiprocessing.Queue): A queue from the `multiprocessing` module where log 
                                       records are placed by other processes. The queue is 
                                       used to communicate log messages to this logger process.
        log_file_path_func (callable): A function that generates the log file path dynamically. 
                                       This function should return the path where the log file 
                                       should be saved.

    Raises:
        Exception: If an error occurs while handling log records (e.g., during logging or 
                   record retrieval from the queue).

    Notes:
        - The log records are formatted and outputted with a timestamp, log level, and message.
        - The function runs in an infinite loop, processing log records until a termination 
          signal (`None`) is received from the queue.
        - Logs are written to both a log file and printed to the console for real-time monitoring.
    """
    # Get the root logger, which is the main logging object in the system
    root = logging.getLogger()

    # Create file and stream handlers to output logs to both a log file and the console
    file_handler = logging.FileHandler(log_file_path)
    stream_handler = logging.StreamHandler()

    # Set up a formatter to specify how logs will appear (with timestamp, log level, and message)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # Add the file and stream handlers to the root logger
    root.addHandler(file_handler)
    root.addHandler(stream_handler)

    # Set the logging level to INFO, which will capture all INFO, WARNING, ERROR, and CRITICAL logs
    root.setLevel(logging.INFO)

    # Start an infinite loop to continually check the queue for log records
    while True:
        try:
            # Retrieve a log record from the queue
            record = queue.get()

            # If the record is None, this is a termination signal, so break the loop
            if record is None:
                break

            # Get the logger for the specific record's name and process the log message
            logger = logging.getLogger(record.name)
            logger.handle(record)

        except Exception:
            # If an error occurs while processing the log record, print the error traceback to stderr
            import sys, traceback
            print('Error in logger process', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)


def setup_subprocess_logger(shared_queue: multiprocessing.Queue) -> logging.Logger:
    """
    Sets up a logger for use in multiprocessing contexts. This logger sends log records 
    to a shared queue, which can then be processed by a separate logging process. The 
    logger is configured with a handler that places log records into the shared queue 
    and sets the logging level to INFO.

    Args:
        shared_queue (multiprocessing.Queue): The queue to which log records will be 
                                               sent. This queue will be used by a separate 
                                               logger process to handle and output logs.

    Returns:
        logging.Logger: The configured logger instance that can be used to log messages 
                        from the current process.

    Notes:
        - This function is typically called by worker processes in a multiprocessing setup.
        - The logger is configured to use a `QueueHandler`, which ensures that log records 
          are passed through the shared queue to a separate logging process.
    """
    # Create a logger specifically for multiprocessing tasks
    logger = logging.getLogger("subprocess_logger")
    
    # Add a handler that sends log records to the shared queue
    logger.addHandler(QueueHandler(shared_queue))
    
    # Set the logging level to INFO, so it captures INFO and higher severity logs
    logger.setLevel(logging.INFO)

    # Return the configured logger to be used in the calling process
    return logger



def get_shared_queue_and_logger_instances() -> Tuple[multiprocessing.Queue, logging.Logger]:
    """
    Creates and returns a shared logging queue and a logger instance for multiprocessing.

    This function sets up a shared `multiprocessing.Queue` for logging purposes, which can 
    be used by multiple processes to send log records to a central logging process. It 
    also creates and configures a logger that uses a `QueueHandler` to send log records 
    to the shared queue. The logger is set to the `INFO` logging level by default.

    Returns:
        Tuple[multiprocessing.Queue, logging.Logger]: A tuple containing:
            - A shared `multiprocessing.Queue` that can be used for logging in multiple processes.
            - A `logging.Logger` instance that uses the shared queue to log messages.

    Example:
        shared_queue, logger = get_shared_queue_and_logger_instances()

    Notes:
        - The function initializes the logger to use the `INFO` level and attaches a 
          `QueueHandler` that directs logs to the shared queue.
        - The shared queue is initialized with a maximum size of `-1`, meaning it can grow dynamically.

    """
    # create shared queue
    shared_queue = multiprocessing.Queue(-1)
    # create a logger
    shared_logger = logging.getLogger("shared_logger")
    # add a handler that uses the shared queue
    shared_logger.addHandler(QueueHandler(shared_queue))
    # set level
    shared_logger.setLevel(logging.INFO)

    return shared_queue, shared_logger


