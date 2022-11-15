class HashNode:
    def __init__(self, key, value, next=None):
        self.key = key
        self.value = value
        self.next = next

    def __str__(self):
        if self.key != None and self.value != None:
            return "key = " + str(self.key) + " and value = " + str(self.value)
        else:
            return None


class Hashtabell:
    def __init__(self, elements):
        self.size = elements * 2 + 3
        self.slots = [None] * self.size
        self.krockantal = 0

    def hashKey(self, key):
        """Tagen från föreläsningsanteckningar"""
        result = 0
        for i in key:
            result = result * 32 + ord(i)
        return result % self.size

    def put(self, key, value):
        index = self.hashKey(key)
        if self.slots[index] == None:
            self.slots[index] = HashNode(key, value)
        else:
            self.krockantal += 1
            krock = self.slots[index]
            self.slots[index] = HashNode(key, value)
            self.slots[index].next = krock

    def search(self, key):
        index = self.hashKey(key)
        currNode = self.slots[index]

        while currNode != None:
            if currNode.key == key:
                return currNode
            currNode = currNode.next

        raise KeyError


def storeArtist(a):
    with open("unique_tracks.txt", "r", encoding="ISO-8859-1") as fil:
        for rad in fil:
            attr = rad.split('<SEP>')
            a.put(attr[0], attr[3])


def main():
    a = Hashtabell(1000000)
    storeArtist(a)
    print(a.search("TRMMJXK128F93315E9"))
    print("Antal krockar:", a.krockantal)


if __name__ == '__main__':
    main()