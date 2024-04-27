import os

class MessageLogger:
    """ Class MessageLogger - to store, print, pop up and write messages including error, warning and notice
    Public:
        <func> MessageLogger - construction method
        <func> collect_error - collect error message and append to error message list
        <func> collect_warning - collect warning message and append to warning message list
        <func> print_error - print error messages in the console
        <func> print_warning - print warning messages in the console
        <func> print_message - print notice message in the console
        <func> popup_error - show in a message box error messages
        <func> popup_warning - show in a message box warning messages
        <func> popup_message - show given message in a message box
        <func> write_error - write into the log file error messages
        <func> write_warning - write into the log file warning messages
        <func> write_message - write given message into the log file
        <func> clear_error - clear all error messages stored in the error list
        <func> clear_warning - clear all warning messages stored in the warning list
    Private:
        <attr> log_file - the path of the log file
        <attr> error_messages - a list to store error messages
        <attr> warning_messages - a list to store warning messages
    """
    
    def __init__(self, log_file="log_file.txt", log_rewrite=True):
        """ construction method 
        ----------
        Parameters
        - log_file - the path of the log file (or the parent folder of the log file)
        - log_rewrite - not yet used
        """
        # initialize by input parameters
        if os.path.isdir(log_file):
            self.__log_file = os.path.join(log_file, "log_file.txt")
        else:
            self.__log_file = log_file
            
        # delete old log_file
        f = open(self.__log_file, "w")
        f.close()
        
        # initialize message lists
        self.__error_messages = []
        self.__warning_messages = []
        
        
    def collect_error(self, error_message):
        """ collect error message and append to error message list
        ----------
        Parameters
        - string error_message - error message to be added to the error list
        """
        self.__error_messages.append(error_message)
        
        
    def collect_warning(self, warning_message):
        """ collect warning message and append to warning message list
        ----------
        Parameters
        - string warning_message - warning message to be added to the warning list
        """
        self.__warning_messages.append(warning_message)
        
        
    def print_error(self, which="all"):
        """ print all error messages if which=="all", otherwise latest error message
        ----------
        Parameters
        - which - print all error messages if which=="all", otherwise latest error
            Default: "all"
        """
        if len(self.__error_messages) == 0:
            return
        
        print("\n===== error message(s) =====")
        if which == "all":  # print all messages
            for msg in self.__error_messages:
                print(f"ERROR: {msg}")
        else:  # print only the latest one
            print(f"ERROR: {self.__error_messages[-1]}")
        print("===== error message(s) =====")
        
        
    def print_warning(self, which="all"):
        """ print all warning messages if which=="all", otherwise latest warning message
        ----------
        Parameters
        - which - print all warning messages if which=="all", otherwise latest warning
            Default: "all"
        """
        if len(self.__warning_messages) == 0:
            return
        
        print("\n===== warning message(s) =====")
        if which == "all":  # print all messages
            for msg in self.__warning_messages:
                print(f"WARNING: {msg}")
        else:  # print only the latest one
            print(f"WARNING: {self.__error_messages[-1]}")
        print("===== warning message(s) =====")

    
    def print_message(self, message):
        """ print notice message in the console 
        ----------
        Parameters
        - message - message to print in the console, a string or a list of string
        """
        if isinstance(message, list):
            print("\n")
            for msg in message:
                print(message)
        else:
            print(message)

        
    def write_error(self, which="all", mode="w"):
        """ write into the log file all error messages if which=="all", otherwise latest error message
        ----------
        Parameters
        - which - write error messages if which=="all", otherwise latest error
            Default: "all"
        - mode - the mode opening the file, "a" for append, "w" for write, and so on
            Default: "w"
        """
        if len(self.__error_messages) == 0:
            return
        
        log_file = open(self.__log_file, mode)
        if which == "all":
            for msg in self.__error_messages:
                log_file.write("ERROR:" + msg + "\n")
        log_file.close()
        
        
    def write_warning(self, which="all", mode="w"):
        """ write into the log file all warning messages if which=="all", otherwise latest warning message
        ----------
        Parameters
        - which - write warning messages if which=="all", otherwise latest warning
            Default: "all"
        - mode - the mode opening the file, "a" for append, "w" for write, and so on
            Default: "w"
        """
        if len(self.__warning_messages) == 0:
            return
        
        log_file = open(self.__log_file, mode)
        if which == "all":
            for msg in self.__warning_messages:
                log_file.write("Warning:" + msg + "\n")
        log_file.close()
        
        
    def write_message(self, message, mode="a"):
        """ write message into the log file 
        ----------
        Parameters
        - message - the message to be written into the log file
        - mode - the mode opening the file, "a" for append, "w" for write, and so on
            Default: "a"
        """
        try:
            log_file = open(self.__log_file, mode)
        except:
            log_file = open(self.__log_file, "w")
        log_file.write(message + "\n")
        log_file.close()
        
        
    def clear_error(self):
        """ clear all error messages stored in the error list """
        self.__error_messages = []
        
        
    def clear_warning(self):
        """ clear all warning messages stored in the warning list """
        self.__warning_messages = []