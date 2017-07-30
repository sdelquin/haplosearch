from .exceptions import HaploException
import re
from .utils import listtoranges

# definición de las bases
ADENINA = "A"
TIMINA = "T"
GUANINA = "G"
CITOSINA = "C"
BLANCO = "-"
# códigos de la IUPAC
ADEoGUA = "R"
CIToTIM = "Y"
GUAoCIT = "S"
ADEoTIM = "W"
GUAoTIM = "K"
ADEoCIT = "M"
CIToGUAoTIM = "B"
ADEoGUAoTIM = "D"
ADEoCIToTIM = "H"
ADEoCIToGUA = "V"
ANYBASE = "N"

# definición de los cambios
TRANSITION = "T"
TRANSVERSION = "V"
INSERTION = "I"
DELETION = "D"
MISSING = "M"
HETEROPLASMY = "H"


def split_mutation_population(mutation):
    """
    Método que separa una mutación en tipo, posición y bases. Se devuelve en un
    diccionario.
    Para el caso de las transiciones y transversiones siempre se devolverá un
    único tipo, posición y base; para el caso de las deleciones e inserciones
    se podrían devolver más de una base, lo que llevaría a que la posición que
    se devuelve es la correspondiente a la primera base.
    Se utiliza la nomenclatura de genética de poblaciones.
    """
    # diccionario que vamos a retornar
    d = {}
    # el elemento es una transición
    m = re.compile(r"^\d+$").match(mutation)
    if (m):
        d["type"] = TRANSITION
        d["position"] = int(m.group(0))
        d["base"] = ""
    else:
        # el elemento es una transversión
        m = re.compile(r"^(\d+)([ATGCRYSWKMBDHVN])$").match(mutation)
        if (m):
            d["type"] = TRANSVERSION
            d["position"] = int(m.group(1))
            d["base"] = m.group(2)
        else:
            # el elemento es una inserción
            m = re.compile(r"^(\d+)i([ATGC]+)$").match(mutation)
            if (m):
                d["type"] = INSERTION
                d["position"] = int(m.group(1))
                d["base"] = m.group(2)
            else:
                # el elemento es una deleción
                m = re.compile(r"^(\d+)d([ATGC]+)$").match(mutation)
                if (m):
                    d["type"] = DELETION
                    d["position"] = int(m.group(1))
                    d["base"] = m.group(2)
                else:
                    d["type"] = ""
                    d["position"] = ""
                    d["base"] = ""
    # retornamos el diccionario
    return d


def split_mutation_forensic(mutation):
    """
    Método que separa una mutación en tipo, posición y bases. Se devuelve en un
    diccionario.
    Para el caso de las transiciones y transversiones siempre se devolverá un
    único tipo, posición y base; para el caso de las deleciones e inserciones
    se podrían devolver más de una base, lo que llevaría a que la posición que
    se devuelve es la correspondiente a la primera base.
    Se utiliza la nomenclatura forense.
    """
    # diccionario que vamos a retornar
    d = {}
    # el elemento puede ser una transición, una transversión ó una
    # heteroplasmia, pero a nuestros efectos, la consideramos como transversión
    m = re.compile(r"^(\d+)([ATGCRYSWKMBDHVN])$").match(mutation)
    if (m):
        d["type"] = TRANSVERSION
        d["position"] = int(m.group(1))
        d["base"] = m.group(2)
    else:
        # el elemento es una inserción
        m = re.compile(r"^(\d+)\.(\d+)([ATGC])$").match(mutation)
        if (m):
            d["type"] = INSERTION
            d["position"] = int(m.group(1)) + int(m.group(2)) - 1
            d["base"] = m.group(3)
        else:
            # el elemento es una deleción
            m = re.compile(r"^(\d+)(del|d)$").match(mutation)
            if (m):
                d["type"] = DELETION
                d["position"] = int(m.group(1))
                # no conocemos la base, así que ponemos missing
                d["base"] = MISSING
            else:
                d["type"] = ""
                d["position"] = ""
                d["base"] = ""
    # retornamos el diccionario
    return d


class Adn:
    """
    Clase que representa una secuencia de ADN.
    """

    def __init__(self, id, sec):
        """
        Constructor de la clase.
        Se construye un objeto de ADN a partir de la secuencia de nucleotidos
        Parámetros:
            id: identificación de la secuencia
            sec: cadena de caracteres con las bases de la secuencia
        """
        self.id = id
        self.sec = list(sec[:])
        self.posref = []
        self.mutations = []

    @classmethod
    def from_haplotype_population(cls, id, ref, haplotipo):
        """
        Constructor "alternativo" de la clase.
        Se construye un objeto de ADN a partir del haplotipo sobre una
        referencia y una posición base.
        Se utiliza la nomenclatura de genética de poblaciones.
        Parámetros:
          id: identificación de la secuencia
          ref: secuencia de referencia (normalmente rCRS)
          haplotipo: haplotipo de la nueva secuencia
        """
        # lista de mutaciones, en principio vacía
        mutations = []
        # en principio, la secuencia es la misma que la referencia
        sec = ref[:]
        # si el haplotipo coincide con la denominación de la referencia es que
        # es la misma y no hay cambios
        if (ref.id != haplotipo):
            # dividimos el haplotipo por los espacios
            haplo_lista = haplotipo.split()
            # recorremos los diferentes elementos del haplotipo
            for h in haplo_lista:
                # separamos la mutación
                m = split_mutation_population(h)
                # guardamos la mutación
                mutations.append(m)
                # el elemento es una transición
                if (m["type"] == TRANSITION):
                    pos = ref.posref.index(int(m["position"]))
                    if (ref[pos] == TIMINA):
                        m["base"] = sec[pos] = CITOSINA
                    elif (ref[pos] == CITOSINA):
                        m["base"] = sec[pos] = TIMINA
                    elif (ref[pos] == ADENINA):
                        m["base"] = sec[pos] = GUANINA
                    elif (ref[pos] == GUANINA):
                        m["base"] = sec[pos] = ADENINA
                # el elemento es una transversión
                elif (m["type"] == TRANSVERSION):
                    pos = ref.posref.index(int(m["position"]))
                    base = m["base"]
                    sec[pos] = base
                # el elemento es una inserción
                elif (m["type"] == INSERTION):
                    pos = ref.posref.index(int(m["position"])) + 1
                    # iteramos por si hay más de una base
                    for b in m["base"]:
                        sec[pos] = b
                        pos = pos + 1
                elif (m["type"] == DELETION):
                    pos = ref.posref.index(int(m["position"]))
                    # iteramos por si hay más de una base
                    for b in m["base"]:
                        sec[pos] = BLANCO
                        pos = pos + 1
                # si entramos por este último else es que lo que aparece en la
                # entrada es incorrecto
                else:
                    raise HaploException(f"Unknown notation: {h}")
        # creamos el nuevo objeto
        s = cls(id, "".join(sec))
        # le asociamos las mutaciones
        s.mutations = mutations
        # retornamos la secuencia
        return s

    @classmethod
    def from_haplotype_forensic(cls, id, ref, haplotipo):
        """
        Constructor "alternativo" de la clase.
        Se construye un objeto de ADN a partir del haplotipo sobre una
        referencia y una posición base.
        Se utiliza la nomenclatura forense.
        Parámetros:
          id: identificación de la secuencia
          ref: secuencia de referencia (normalmente rCRS)
          haplotipo: haplotipo de la nueva secuencia
        """
        # lista de mutaciones, en principio vacía
        mutations = []
        # en principio, la secuencia es la misma que la referencia
        sec = ref[:]
        # si el haplotipo coincide con la denominación de la referencia es que
        # es la misma y no hay cambios
        if (ref.id != haplotipo):
            # dividimos el haplotipo por los espacios
            haplo_lista = haplotipo.split()
            # recorremos los diferentes elementos del haplotipo
            for h in haplo_lista:
                # el elemento lo tratamos como transversión
                m = re.compile(r"^(\d+)([ATGCRYSWKMBDHVN])$").match(h)
                if (m):
                    pos = ref.posref.index(int(m.group(1)))
                    base = m.group(2)
                    sec[pos] = base
                else:
                    # el elemento es una inserción
                    m = re.compile(r"^(\d+)\.(\d+)([ATGC])$").match(h)
                    if (m):
                        base = m.group(3)
                        pos = ref.posref.index(int(m.group(1))) + \
                            int(m.group(2))
                        sec[pos] = base
                    else:
                        # el elemento es una deleción
                        m = re.compile(r"^(\d+)(del|d)$").match(h)
                        if (m):
                            pos = ref.posref.index(int(m.group(1)))
                            sec[pos] = BLANCO
                        # si entramos por este último else es que lo que
                        # aparece en la entrada es incorrecto
                        else:
                            raise HaploException(f"Unknown notation: {h}")
        # creamos el nuevo objeto
        s = cls(id, "".join(sec))
        # le asociamos las mutaciones
        s.mutations = mutations
        # retornamos la secuencia
        return s

    def build_positions(self, posbase):
        """
        Función para construir las posiciones en base a una secuencia de
        referencia y una posición inicial.
        Parámetros:
          ref: secuencia de referencia
          posbase: posición base de la secuencia de referencia
        Retorna:
          i: lista de posiciones en base a la secuencia
        """
        # inicializar la lista a vacío
        self.posref = []
        # contador auxiliar
        c = posbase
        # colocamos la posición inicial
        self.posref.append(c)
        # recorremos el resto de elementos
        for s in self.sec[1:]:
            if (s == BLANCO):
                self.posref.append(c)
            else:
                c = c + 1
                self.posref.append(c)

    # Imprimir la secuencia
    def __str__(self):
        return ">%s\n%s" % (self.id, "".join(self.sec))

    # Imprimir la secuencia sin huecos
    def versinhuecos(self):
        return ">%s\n%s" % (self.id, "".join(self.sec).replace(BLANCO, ""))

    def getseq_asstring(self):
        """
        Función que devuelve la secuencia de ADN como un string.
        """
        return "".join(self.sec)

    def __getitem__(self, position):
        """
        Obtener una base de la secuencia.
        Parámetros:
          position: posición que se quiere buscar
        """
        return self.sec[position]

    def __len__(self):
        """
        Obtener la longitud de la secuencia.
        """
        return len(self.sec)

    def haplotype_population(self, ref, haps3digits):
        """
        Obtener el haplotipo de la secuencia actual en comparación con una de
        referencia.
        Se utiliza la nomenclatura de genética de poblaciones.
        Parámetros:
            ref: secuencia de referencia.
            haps3digits: indica si los haplotipos se muestran con 3 dígitos o
            normales.
        Retorna:
            Una cadena con la representación del haplotipo.
        """
        # posición de la secuencia
        i = 0
        # lista para guardar los valores missing
        missing = []
        # hay que llevar cuenta del último cambio registrado
        ultimo_cambio = -1
        # haplotipo
        h = ""
        # recorremos los elementos de las dos secuencias
        # se presupone que tienen el mismo tamaño
        while (i < len(self)):
            # posición actual de referencia
            pos = ref.posref[i]
            # TRANSICIONES
            if (((ref[i] == CITOSINA) and (self[i] == TIMINA)) or
                    ((ref[i] == TIMINA) and (self[i] == CITOSINA)) or
                    ((ref[i] == ADENINA) and (self[i] == GUANINA)) or
                    ((ref[i] == GUANINA) and (self[i] == ADENINA))):
                if (haps3digits):
                    h = h + (" %03d" % (pos))
                else:
                    h = h + (" %d" % (pos))
                ultimo_cambio = TRANSITION
            # TRANSVERSIONES
            elif (((ref[i] == ADENINA) and (self[i] == TIMINA)) or
                  ((ref[i] == TIMINA) and (self[i] == ADENINA)) or
                  ((ref[i] == ADENINA) and (self[i] == CITOSINA)) or
                  ((ref[i] == CITOSINA) and (self[i] == ADENINA)) or
                  ((ref[i] == GUANINA) and (self[i] == TIMINA)) or
                  ((ref[i] == TIMINA) and (self[i] == GUANINA)) or
                  ((ref[i] == GUANINA) and (self[i] == CITOSINA)) or
                  ((ref[i] == CITOSINA) and (self[i] == GUANINA))):
                if (haps3digits):
                    h = h + (" %03d%s" % (pos, self[i]))
                else:
                    h = h + (" %d%s" % (pos, self[i]))
                ultimo_cambio = TRANSVERSION
            # HETEROPLASMIAS
            elif ((self[i] == ADEoGUA) or
                  (self[i] == ADEoGUA) or
                  (self[i] == CIToTIM) or
                  (self[i] == GUAoCIT) or
                  (self[i] == ADEoTIM) or
                  (self[i] == GUAoTIM) or
                  (self[i] == ADEoCIT) or
                  (self[i] == CIToGUAoTIM) or
                  (self[i] == ADEoGUAoTIM) or
                  (self[i] == ADEoCIToTIM) or
                  (self[i] == ADEoCIToGUA) or
                  (self[i] == ANYBASE)):
                if (haps3digits):
                    h = h + (" %03d%s" % (pos, self[i]))
                else:
                    h = h + (" %d%s" % (pos, self[i]))
                ultimo_cambio = HETEROPLASMY
            # INSERCIONES
            elif (((ref[i] == BLANCO) and (self[i] == TIMINA)) or
                  ((ref[i] == BLANCO) and (self[i] == ADENINA)) or
                  ((ref[i] == BLANCO) and (self[i] == GUANINA)) or
                  ((ref[i] == BLANCO) and (self[i] == CITOSINA))):
                # si el último cambio no fue una inserción, empezamos una nueva
                if (ultimo_cambio != INSERTION):
                    if (haps3digits):
                        h = h + (" %03di%s" % (pos, self[i]))
                    else:
                        h = h + (" %di%s" % (pos, self[i]))
                # si el último cambio fue una inserción, la agrupamos con la
                # anterior
                else:
                    h = h + ("%s" % (self[i]))
                ultimo_cambio = INSERTION
            # DELECIONES
            elif (((ref[i] == TIMINA) and (self[i] == BLANCO)) or
                  ((ref[i] == ADENINA) and (self[i] == BLANCO)) or
                  ((ref[i] == GUANINA) and (self[i] == BLANCO)) or
                  ((ref[i] == CITOSINA) and (self[i] == BLANCO))):
                # si el último cambio no fue una deleción, empezamos una nueva
                if (ultimo_cambio != DELETION):
                    if (haps3digits):
                        h = h + (" %03dd%s" % (pos, ref[i]))
                    else:
                        h = h + (" %dd%s" % (pos, ref[i]))
                else:
                    h = h + ("%s" % (ref[i]))
                ultimo_cambio = DELETION
            # MISSING
            elif ((self[i] != TIMINA) and
                  (self[i] != ADENINA) and
                  (self[i] != GUANINA) and
                  (self[i] != CITOSINA) and
                  (self[i] != BLANCO)):
                missing.append(pos)
                ultimo_cambio = MISSING
            else:
                ultimo_cambio = -1

            # incrementar el contador
            i = i + 1

        # añadir los valores missing en el caso de que no sea vacío
        if (len(missing) > 0):
            h = h + "  MISSING: " + listtoranges(missing)
        # si el haplotipo es vacío coincide con el CRS
        if (h == ""):
            h = ref.id
        # retornamos el haplotipo sin espacios antes y después
        return h.strip()

    def haplotype_forensic(self, ref, haps3digits, deletions_as_d):
        """
        Obtener el haplotipo de la secuencia actual en comparación con una de
        referencia.
        Se utiliza la nomenclatura forense.
        Parámetros:
          ref: secuencia de referencia.
            haps3digits: indica si los haplotipos se muestran con 3 dígitos o
                normales.
            deletions_as_d: controla si las deleciones se muestran como "d" ó
                "del" sólo en nomenclatura forense
        Retorna:
          Una cadena con la representación del haplotipo.
        """
        # posición de la secuencia
        i = 0
        # contador para las inserciones múltiples
        ins_counter = 0
        # lista para guardar los valores missing
        missing = []
        # haplotipo
        h = ""
        # cadena de representación de deleciones
        if (deletions_as_d):
            deletions_str = "d"
        else:
            deletions_str = "del"
        # recorremos los elementos de las dos secuencias
        # se presupone que tienen el mismo tamaño
        while (i < len(self)):
            # posición actual de referencia
            pos = ref.posref[i]
            # TRANSICIONES
            if (
                ((ref[i] == CITOSINA) and (self[i] == TIMINA)) or
                ((ref[i] == TIMINA) and (self[i] == CITOSINA)) or
                ((ref[i] == ADENINA) and (self[i] == GUANINA)) or
                ((ref[i] == GUANINA) and (self[i] == ADENINA)) or
                # TRANSVERSIONES
                ((ref[i] == ADENINA) and (self[i] == TIMINA)) or
                ((ref[i] == TIMINA) and (self[i] == ADENINA)) or
                ((ref[i] == ADENINA) and (self[i] == CITOSINA)) or
                ((ref[i] == CITOSINA) and (self[i] == ADENINA)) or
                ((ref[i] == GUANINA) and (self[i] == TIMINA)) or
                ((ref[i] == TIMINA) and (self[i] == GUANINA)) or
                ((ref[i] == GUANINA) and (self[i] == CITOSINA)) or
                ((ref[i] == CITOSINA) and (self[i] == GUANINA)) or
                # HETEROPLASMIAS
                (self[i] == ADEoGUA) or
                (self[i] == ADEoGUA) or
                (self[i] == CIToTIM) or
                (self[i] == GUAoCIT) or
                (self[i] == ADEoTIM) or
                (self[i] == GUAoTIM) or
                (self[i] == ADEoCIT) or
                (self[i] == CIToGUAoTIM) or
                (self[i] == ADEoGUAoTIM) or
                (self[i] == ADEoCIToTIM) or
                (self[i] == ADEoCIToGUA) or
                (self[i] == ANYBASE)
            ):
                if (haps3digits):
                    h = h + (" %03d%s" % (pos, self[i]))
                else:
                    h = h + (" %d%s" % (pos, self[i]))
                ins_counter = 0
            # INSERCIONES
            elif (
                ((ref[i] == BLANCO) and (self[i] == TIMINA)) or
                ((ref[i] == BLANCO) and (self[i] == ADENINA)) or
                ((ref[i] == BLANCO) and (self[i] == GUANINA)) or
                ((ref[i] == BLANCO) and (self[i] == CITOSINA))
            ):
                ins_counter += 1
                if (haps3digits):
                    h = h + (" %03d.%d%s" % (pos, ins_counter, self[i]))
                else:
                    h = h + (" %d.%d%s" % (pos, ins_counter, self[i]))
            # DELECIONES
            elif (
                ((ref[i] == TIMINA) and (self[i] == BLANCO)) or
                ((ref[i] == ADENINA) and (self[i] == BLANCO)) or
                ((ref[i] == GUANINA) and (self[i] == BLANCO)) or
                ((ref[i] == CITOSINA) and (self[i] == BLANCO))
            ):
                if (haps3digits):
                    h = h + (" %03d%s" % (pos, deletions_str))
                else:
                    h = h + (" %d%s" % (pos, deletions_str))
                ins_counter = 0
            # MISSING
            elif (
                (self[i] != TIMINA) and
                (self[i] != ADENINA) and
                (self[i] != GUANINA) and
                (self[i] != CITOSINA) and
                (self[i] != BLANCO)
            ):
                missing.append(pos)
                ins_counter = 0
            else:
                ins_counter = 0

            # incrementar el contador
            i = i + 1

        # añadir los valores missing en el caso de que no sea vacío
        if (len(missing) > 0):
            h = h + "  MISSING: " + listtoranges(missing)
        # si el haplotipo es vacío coincide con el CRS
        if (h == ""):
            h = ref.id
        # retornamos el haplotipo sin espacios antes y después
        return h.strip()

    def align_sequence_population(self, begin_pos, nom_fich):
        """
        Función que alinea una secuencia. El procedimiento es buscar en el
        fichero de entrada todas las inserciones que existan e ir añadiendo
        huecos a la secuencia de entrada, teniendo en cuenta que hay que
        recalcular las posiciones en base a la posición que se pasa.
        Se utiliza la nomenclatura de genética de poblaciones.
        Parámetros:
          begin_pos: posición de comienzo.
          nom_fich: nombre del fichero de entrada.
        """
        # en primer lugar calculamos las posiciones
        self.build_positions(begin_pos)
        # apertura del fichero
        f = open(nom_fich)
        # recorrer el fichero
        for linea in f.readlines():
            # recorremos cada elemento de la línea
            for e in linea.strip().split():
                # buscamos inserciones
                m = re.compile(r"^(\d+)i([ATGC]+)$").match(e)
                if (m):
                    # se suma uno porque se supone que el hueco está en el
                    # siguiente nucleótido
                    pos = self.posref.index(int(m.group(1))) + 1
                    # insertamos tantos huecos como número de bases hayan
                    # especificado en la inserción
                    for i in range(len(m.group(2))):
                        # sólo insertamos el hueco si ahí ya no había otro
                        if (self.sec[pos] != BLANCO):
                            self.sec.insert(pos, BLANCO)
                        pos += 1
                    # hay que recalcular las posiciones
                    self.build_positions(begin_pos)
        # cerramos el fichero
        f.close()

    def align_sequence_forensic(self, begin_pos, nom_fich):
        """
        Función que alinea una secuencia. El procedimiento es buscar en el
        fichero de entrada todas las inserciones que existan e ir añadiendo
        huecos a la secuencia de entrada, teniendo en cuenta que hay que
        recalcular las posiciones en base a la posición que se pasa.
        Se utiliza la nomenclatura forense.
        Parámetros:
          begin_pos: posición de comienzo.
          nom_fich: nombre del fichero de entrada.
        """
        # en primer lugar calculamos las posiciones
        self.build_positions(begin_pos)
        # apertura del fichero
        f = open(nom_fich)
        # recorrer el fichero
        for linea in f.readlines():
            # recorremos cada elemento de la línea
            for e in linea.strip().split():
                # buscamos inserciones
                m = re.compile(r"^(\d+)\.(\d+)([ATGC])$").match(e)
                if (m):
                    # si el desplazamiento es mayor que 0 significa que estamos
                    # una inserción múltiple, en otro caso es la primera
                    # inserción lo que conlleva sumar uno
                    pos = self.posref.index(int(m.group(1))) + int(m.group(2))
                    # sólo insertamos el hueco si ahí ya no había otro
                    if (self.sec[pos] != BLANCO):
                        self.sec.insert(pos, BLANCO)
                    # hay que recalcular las posiciones
                    self.build_positions(begin_pos)
        # cerramos el fichero
        f.close()


def sec2hap(nom_fichero_entrada, nom_fichero_salida, nomenclature):
    """
    Función que lee secuencias de un fichero de entrada y construye
    los haplotipos asociados en base a una referencia, y además los
    escribe en un fichero de salida.
    Parámetros:
    nom_fichero_entrada: nombre del fichero de datos de entrada (haplotipos)
    nom_fichero_salida: nombre del fichero de datos de salida (secuencias)
    La estructura del fichero de secuencias es la siguiente:
    START: x  (indica la posición en el que se empieza a contar la secuencia
    de referencia)
    >[IDENTIFICACION DE LA SECUENCIA DE REFERENCIA]
    [SECUENCIA DE REFERENCIA]
    >[IDENTIFICACIÓN DE LA PRIMERA SECUENCIA]
    [PRIMERA SECUENCIA]
    >[IDENTIFICACIÓN DE LA SEGUNDA SECUENCIA]
    [SEGUNDA SECUENCIA]
    ...
    nomenclature: tipo de nomenclatura (poblaciones ó forense)
    """
    # apertura de los ficheros
    fichero_entrada = open(nom_fichero_entrada, "r")
    fichero_salida = open(nom_fichero_salida, "w")
    # en la primera fila se encuentra la posición base
    l = fichero_entrada.readline().strip()
    g = re.compile(r"^START *: *(\d+) *(\*|d|\*d|d\*)? *$").match(l)
    # variable que controla si vamos a sacar las posiciones de los haplotipos
    # con 3 dígitos o normales
    haps3digits = False
    # variable que controla si las deleciones se muestran como "d" ó "del" sólo
    # en nomenclatura forense
    deletions_as_d = False
    if (g):
        posbase = int(g.group(1))
        options = g.group(2)
        if (options):
            if ("*" in options):
                haps3digits = True
            if ("d" in options):
                deletions_as_d = True
    else:
        raise HaploException(
            "Base position is not defined on input file"
        )
    # inicializamos el número de línea
    num_linea = 2
    # en la primera línea está la identificación de la secuencia de referencia
    m = re.compile(r"^ *> *(.+)$").match(fichero_entrada.readline().rstrip())
    # controlamos un posible error de sintaxis
    if (not m):
        raise HaploException(f"Syntax error on input file [Line: {num_linea}]")
    id = m.group(1)
    # en la segunda línea está la secuencia de referencia en sí misma
    num_linea = 3
    m = re.compile(r"^ *([AGTC-]+) *$").match(
        fichero_entrada.readline().strip().upper()
    )
    # controlamos un posible error de sintaxis
    if (not m):
        raise HaploException(f"Syntax error on input file [Line: {num_linea}]")
    sec = m.group(1)
    # construimos la secuencia de referencia
    ref = Adn(id, sec)
    # construimos las posiciones de referencia
    ref.build_positions(posbase)
    # escribimos la posición base en el fichero de salida
    fichero_salida.write("START: %d\n" % (posbase))
    # escribimos la secuencia de referencia en el fichero de salida
    fichero_salida.write("%s\n" % (ref.versinhuecos()))
    # procesamiento del resto del fichero
    num_linea = 4
    # recorremos el fichero de entrada
    while (True):
        # los datos vienen emparejados
        # la primera línea es la identificación de la secuencia
        linea = fichero_entrada.readline().strip()
        if (linea == ""):
            break
        m = re.compile(r"^ *> *(.+)$").match(linea)
        # controlamos un posible error de sintaxis
        if (not m):
            raise HaploException(
                f"Syntax error on input file [Line: {num_linea}]"
            )
        # si todo ha ido bien tenemos el identificador en el primer grupo
        id = m.group(1)
        # la segunda línea es la secuencia de adn
        num_linea = num_linea + 1
        m = re.compile(r"^ *([AGTCRYSWKMBDHVN-]+) *$").match(
            fichero_entrada.readline().strip().upper()
        )
        # controlamos un posible error de sintaxis
        if (not m):
            raise HaploException(
                f"Syntax error on input file [Line: {num_linea}]"
            )
        # si todo ha ido bien tenemos la secuencia en el primer grupo
        sec = m.group(1)
        # controlamos que la secuencia de entrada tenga el mismo número de
        # bases que la secuencia de referencia
        if (len(sec) != len(ref)):
            raise HaploException(
                f"Sequence has not the same number of bases as the reference "
                "sequence [Line: {num_linea}]"
            )
        # construimos una nueva secuencia
        a = Adn(id, sec)
        # escribimos la secuencia en el fichero de salida
        try:
            if (nomenclature == "POP"):
                h = a.haplotype_population(ref, haps3digits)
            else:
                h = a.haplotype_forensic(ref, haps3digits, deletions_as_d)
            fichero_salida.write(">%s\n%s\n" % (a.id, h))
        except HaploException as err:
            msg = f"Data error on input file [Line: {num_linea}]"
            raise HaploException(msg + "\n" + err.args[0])
        # incrementamos el número de línea
        num_linea = num_linea + 1
    # cerramos los ficheros
    fichero_entrada.close()
    fichero_salida.close()


def hap2sec(nom_fichero_entrada, nom_fichero_salida, nomenclature):
    """
    Función que lee haplotipos de un fichero de entrada y construye
    las secuencias asociadas en base a una referencia, y además las
    escribe en un fichero de salida.
    Parámetros:
    nom_fichero_entrada: nombre del fichero de datos de entrada (haplotipos)
    nom_fichero_salida: nombre del fichero de datos de salida (secuencias)

    La estructura del fichero de haplotipos es la siguiente:
    START: x  (indica la posición en el que se empieza a contar la secuencia de
    referencia)
    >[IDENTIFICACION DE LA SECUENCIA DE REFERENCIA]
    [SECUENCIA DE REFERENCIA]
    >[IDENTIFICACIÓN DE LA PRIMERA SECUENCIA]
    [HAPLOTIPO DE LA PRIMERA SECUENCIA]
    >[IDENTIFICACIÓN DE LA SEGUNDA SECUENCIA]
    [HAPLOTIPO DE LA SEGUNDA SECUENCIA]
    ...
    nomenclature: tipo de nomenclatura (poblaciones ó forense)
    """
    # apertura de los ficheros
    fichero_entrada = open(nom_fichero_entrada, "r")
    fichero_salida = open(nom_fichero_salida, "w")
    # en la primera fila se encuentra la posición base
    l = fichero_entrada.readline().strip()
    g = re.compile(r"^START: *(\d+)").match(l)
    if (g):
        posbase = int(g.group(1))
    else:
        raise HaploException("Base position is not defined on input file")
    # inicializamos el número de línea
    num_linea = 2
    # en la segunda línea está la identificación de la secuencia de referencia
    m = re.compile(r"^ *> *(.+)$").match(fichero_entrada.readline().rstrip())
    # controlamos un posible error de sintaxis
    if (not m):
        raise HaploException(
            f"Syntax error on input file [Line: {num_linea}]"
        )
    id = m.group(1)
    # en la tercera línea está la secuencia de referencia en sí misma
    num_linea = 3
    m = re.compile(r"^ *([AGTC]+) *$").match(
        fichero_entrada.readline().strip().upper()
    )
    # controlamos un posible error de sintaxis
    if (not m):
        raise HaploException(
            f"Syntax error on input file [Line: {num_linea}]"
        )
    sec = m.group(1)
    # construimos la secuencia de referencia
    ref = Adn(id, sec)
    # cerrar temporalmente el fichero de entrada
    fichero_entrada.close()
    # alinear la secuencia de referencia
    if (nomenclature == "POP"):
        ref.align_sequence_population(posbase, nom_fichero_entrada)
    else:
        ref.align_sequence_forensic(posbase, nom_fichero_entrada)
    # abrimos de nuevo el fichero de entrada e ignoramos las tres primeras
    # líneas
    fichero_entrada = open(nom_fichero_entrada, "r")
    fichero_entrada.readline()
    fichero_entrada.readline()
    fichero_entrada.readline()
    # escribimos la posición base en el fichero de salida
    fichero_salida.write("START: %d\n" % (posbase))
    # escribimos la secuencia de referencia en el fichero de salida
    fichero_salida.write("%s\n" % (ref))
    # procesamiento del resto del fichero
    num_linea = 4
    while (True):
        # los datos vienen emparejados
        # la primera línea es la identificación de la secuencia
        linea = fichero_entrada.readline().strip()
        if (linea == ""):
            break
        m = re.compile(r"^ *> *(.+)$").match(linea)
        # controlamos un posible error de sintaxis
        if (not m):
            raise HaploException(
                f"Syntax error on input file [Line: {num_linea}]"
            )
        # si todo ha ido bien tenemos el identificador en el primer grupo
        id = m.group(1)
        num_linea = num_linea + 1
        # la segunda línea es el haplotipo
        hap = fichero_entrada.readline().strip()
        # construimos una nueva secuencia
        try:
            if (nomenclature == "POP"):
                a = Adn.from_haplotype_population(id, ref, hap)
            else:
                a = Adn.from_haplotype_forensic(id, ref, hap)
        except HaploException as err:
            msg = f"Data error on input file [Line: {num_linea}]"
            raise HaploException(msg + "\n" + err.args[0])
        # escribimos la secuencia en el fichero de salida
        fichero_salida.write("%s\n" % (a))
        # incrementamos el número de línea
        num_linea = num_linea + 1
    # cerramos los ficheros
    fichero_entrada.close()
    fichero_salida.close()


def manage_haplosearch(
    inputfile_path, outputfile_path, nomenclature, operation
):
    if operation == "S2H":
        sec2hap(inputfile_path, outputfile_path, nomenclature)
    else:
        hap2sec(inputfile_path, outputfile_path, nomenclature)
