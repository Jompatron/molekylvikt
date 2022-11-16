from molgrafik import *
from hashtest import *

# Syntax:
#<formel>::= <mol> \n
#<mol>   ::= <group> | <group><mol>
#<group> ::= <atom> |<atom><num> | (<mol>) <num>
#<atom>  ::= <LETTER> | <LETTER><letter>
#<LETTER>::= A | B | C | ... | Z
#<letter>::= a | b | c | ... | z
#<num>   ::= 2 | 3 | 4 | ...

period = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
        "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
        "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
        "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        "Fl", "Lv"]
small = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
         'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
big = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
       'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
number = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']


hashadeAtomer = lagraHashtabell(skapaAtomlista())


class LinkedQ:

    class Node:
        #Bestämmer vart de olika pekarna pekar
        def __init__(self, element, next):
            self._next = next
            self.element = element

    #Bestämmer ordningen på noderna
    def __init__(self):
        self._first = None
        self._last = None
        self.size = 0

    def __str__(self):
        return self._first.element

    #Returnerar storleken på den länkade listan
    def __len__(self):
        return self.size

    #Tar bort det element som är först i listan
    def dequeue(self):
        result = self._first.element
        self._first = self._first._next #ändrar pekaren till att nästa element blir det första
        self.size -= 1
        if self.isEmpty():
            self._last = None
        return result

    # Lägger till det element som varit först i listan sist istället
    def enqueue(self, element):
        new = self.Node(element, None)  # Skapar ny nod med element pekar på None
        if self.isEmpty():
            self._first = new  # om kön är tom blir nya noden första noden
        else:
            self._last._next = new  # om det finns i kön läggs noden till sist
        self._last = new
        self.size += 1

    #Kollar om kön är tom
    def isEmpty(self):
        return self.size == 0

    def peek(self):
        if not self.isEmpty():
            return self._first.element
        else:
            return None


class Stack:
    def __init__(self):
        self.stack = []
        self.size = 0

    def store(self, other):
        self.stack.append(other)
        self.size = self.size + 1

    def remove(self):
        last = self.stack.pop()
        self.size = self.size - 1
        return last


class Grammatikfel(Exception):
    pass

class Ruta:
    def __init__(self, atom = "()", num = 1):
        self.atom = atom
        self.num = num
        self.next=None
        self.down=None


def readMolekyl(q, z):
    mol = readGroup(q, z)
    if q.peek() == ".":
        return mol
    elif q.peek() == ")":
        if z.size < 1:
            raise Grammatikfel("Felaktig gruppstart vid radslutet")
        else:
            return mol
    else:
        mol.next = readMolekyl(q, z) 
        
    return mol


def readGroup(q, z):
    rutan = Ruta()
    if q.peek() == ".":
        raise Grammatikfel("Felaktig gruppstart vid radslutet")
    elif q.peek() in number:
        raise Grammatikfel("Felaktig gruppstart vid radslutet")
    elif q.peek().isalpha():
        rutan.atom = readAtom(q)
        if q.peek() == ".":
            return rutan
        if q.peek() in number:
            rutan.num = int(readNum(q))
        

    elif q.peek() == "(":
        z.store(q.dequeue())
        rutan.down = readMolekyl(q, z)
        if q.peek() != ")":
            raise Grammatikfel("Saknad högerparentes vid radslutet")
        if q.peek() == ".":
            raise Grammatikfel("Saknad siffra vid radslutet")
        else:
            z.remove()
            q.dequeue()
            if q.peek() == ".":
                raise Grammatikfel("Saknad siffra vid radslutet")
            rutan.num = int(readNum(q))
    else:
        raise Grammatikfel("Felaktig gruppstart vid radslutet")

    return rutan


def readAtom(q):
    if q.peek() in period:
        first = q.dequeue()
        if q.peek() in small:
            if first + q.peek() in period:
                last = readLetter(q)
                return first + last
            else:
                q.dequeue()
                raise Grammatikfel("Okänd atom vid radslutet ")
        else:
            return first
            
    elif q.peek() in small:
        raise Grammatikfel("Saknad stor bokstav vid radslutet")
    else:
        first = q.dequeue()
        if first + q.peek() in period:
            last = readLetter(q)
            return first + last
        elif first in small:
            raise Grammatikfel("Saknad stor bokstav vid radslutet")
        elif q.peek() in small:
            q.dequeue()
            raise Grammatikfel("Okänd atom vid radslutet")
        else:
            raise Grammatikfel("Okänd atom vid radslutet")


def readLetter(q):
    word = q.peek()
    if word in small:
        q.dequeue()
        return word
    raise Grammatikfel("Fel: Ska följa med liten bokstav: ")


def readNum(q):
    if q.peek() == "0":
        q.dequeue()
        raise Grammatikfel("För litet tal vid radslutet")
    elif q.peek() == "1":
        first = q.dequeue()
        second = q.peek()
        if second in number:
            num = ""
            while q.peek() in number:
                num = num + q.dequeue()
            return num
        else:
            raise Grammatikfel("För litet tal vid radslutet")
    elif q.peek() in number:
        nu = ""
        while q.peek() in number:
            nu = nu + q.dequeue()
        return nu
    else:
        raise Grammatikfel("Saknad siffra vid radslutet")


def storeMolekyl(molekyl):
    q = LinkedQ()
    for letter in molekyl:
        q.enqueue(letter)
    return q


def printRest(q):
    rest = " "
    while not q.isEmpty():
        sak = q.dequeue()
        if sak != ".":
            rest = rest + sak
    return rest


# Inspiration hittad på Internet
# Arvid Larzzon/tilda-labbar på GitHub
def weigh(ruta):
    if ruta.atom == "()":
        if ruta.down is not None:
            x = float(weigh(ruta.down)) * float(ruta.num)

    else:
        x = float(hashadeAtomer.search(ruta.atom).value.vikt) * float(ruta.num)

    if ruta.next is not None:
        y = weigh(ruta.next)
        return x + y
    else:
        return x


def kollaGrammatiken(molekyl):
    molekyl = molekyl + "."
    q = storeMolekyl(molekyl)
    z = Stack()
    try:
        mol = readMolekyl(q, z)
        print("Formeln är syntaktiskt korrekt! Vikten är", weigh(mol))
        return mol
    except Grammatikfel as fel:
        print(str(fel) + printRest(q)) 


def main():
    q = LinkedQ()
    molekyl = input()
    mg = Molgrafik()
    while molekyl != "#":
        mol = kollaGrammatiken(molekyl)
        if mol is str():
            pass
        else:
            #print(mol)
            mg.show(mol)
            molekyl = input()


if __name__ == '__main__':
    main()

