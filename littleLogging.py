# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 19:02:42 2019
@author: solis
Módulo para almacenar mensajes durante la ejecución de un script
Los mensajes se dan de alta con la función append. Los mensajes se graban a
    un fichero de texto con la función dump.
module variables
    __messages. Lista destrings donde se almacenan los mensajes
    __fileNameWihoutExtension. Nombre de fichero de texto donde se grabarán
        los mensajes
    __FILE_EXTENSION. Extensión del fichero log
    __nwrites. Controla el número de veces que se ejecuta dump con éxito
    max_rows. Número máximo de elementos en __messages. Si al hacer un append
        tiene 100 elementos, el programa hace un dump y borra los mensajes
        grabados (control para logs potencialmente de muchos elementos)
"""
from datetime import datetime

__messages = []
__fileNameWihoutExtension = 'app'
__FILE_EXTENSION = '.log'
__nwrites = 0
max_rows = 100


def get_as_list() -> list:
    """
    devuelve el valor de __messages
    """
    return __messages


def get_as_str() -> str:
    """
    devuelve el valor de __messages as str
    """
    return '\n'.join(__messages)


def file_name_get():
	"""
	devuelve el nombre del fichero log
	"""
	return __fileNameWihoutExtension + __FILE_EXTENSION


def append(message: str, toScreen: bool=True):
    """
    append message in __messages
    """
    global __messages
    now = datetime.now()
    if not __messages:
        date = now.strftime('%Y-%m-%d')
        __messages.append(f'{date}')

    if len(__messages) == max_rows:
        dump()
        __messages = []

    time = now.strftime('%HH:%MM:%SS')
    __messages.append(f'{time}: {message}')
    if toScreen:
        print(message)


def dump(fileName: str=__fileNameWihoutExtension, mode: str='w'):
    """
    graba __messages en el fichero fileName
    """
    global __nwrites

    if fileName is not None:
        __fileName = fileName

    if __nwrites == 0:
        mode = 'w'
    else:
        mode = 'a'

    with open(f'{__fileName}{__FILE_EXTENSION}', mode) as f:
        for message in __messages:
            f.write(f'{message}\n')
    __nwrites += 1
