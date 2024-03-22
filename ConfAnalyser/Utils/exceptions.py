class InvalidSourceDataException(Exception):
    """
    Exception to be thrown out in case of invalid data located
    within the source file.
    """
    def __init__(self, message):
        self.message = message

