def handle_uploaded_file(f, path):
    with open(path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    return path


def listtoranges(l):
    """
    Función que convierte una lista de números en sus rangos.
    Parámetros:
        l: lista de números.
    Retorna:
        una cadena con los rangos ó valores sueltos.
    EJEMPLO:
        entrada: [1, 3, 4, 5, 7, 9, 10]
        salida: '1, 3-5, 7, 9-10'
    """
    ini = 0
    fin = 0
    rango = '['
    while (fin < len(l)):
        # mientras el siguiente sea igual al anterior + 1
        while ((fin + 1 < len(l)) and (l[fin + 1] == l[fin] + 1)):
            fin = fin + 1
        # aquí tenemos un rango
        if (ini != fin):
            rango = rango + ('%03d-%03d, ' % (l[ini], l[fin]))
        # aquí tenemos un valor aislado
        else:
            rango = rango + ('%03d, ' % (l[ini]))
        # incrementar final
        fin = fin + 1
        # avanzar inicio
        ini = fin
    # eliminar la última coma y añadir un corchete
    rango = rango[:len(rango) - 2] + ']'
    # retornar el rango
    return rango
