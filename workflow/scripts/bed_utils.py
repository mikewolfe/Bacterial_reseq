class BedEntry(object):
    def __init__(self, field_data=None, names=None, dtypes=None):
        if field_data:
            self.field_data=field_data[:]
        else:
            self.field_data=[]
        if names:
            self.names = {}
            for i, name in enumerate(names):
                self.names[name] = i
        else:
            self.names={"chrm":0, 
                    "start":1, "end":2, "name":3,"score":4, "strand":5,
                    "thickStart":6, "thickEnd":7, "itemRgb":8,
                    "blockCount":9, "blockSizes":10, "blockStarts":11}
        if dtypes:
            self.dtypes=dtypes[:]
        else:
            self.dtypes=[str, int, int, str, float, str,
                         float, float, str, float, str, str]

    def __iter__(self):
        for field in self.field_data:
            yield field

    def __len__(self):
        return len(self.field_data)

    def __repr__(self):
        string = self.field_data[0]
        for field in self.field_data[1:]:
            string += "\t%s"%(field)
        string += "\n"
        return string

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.field_data[key]
        elif isinstance(key, str):
            try:
                name_idx = self.names[key]
            except KeyError:
                raise KeyError("Name does not exist in entry %s"%(key))
            try:
                val = self.field_data[name_idx]
            except IndexError:
                val = "."
            return val
        else:
            raise KeyError("Invalid key %s"%(key))


    def __setitem__(self, key, item):
        if isinstance(key, int):
            if key >= len(self):
                while len(self) < key +1:
                    self.field_data.append(".")
            self.field_data[key] = item
        elif isinstance(key, str):
            try:
                name_idx = self.names[key]
            except KeyError:
                raise KeyError("Name does not exist in entry %s"%(key))
            self.__setitem__(name_idx, item)
        else:
            raise KeyError("Invalid key %s"%(key))

    def parse_from_line(self, line):
        data = []
        for i, (convert, field) in enumerate(zip(self.dtypes, line.rstrip().split("\t"))):
            try:
                datum = convert(field)
            except ValueError:
                datum = str(field)
                self.dtypes[i] = str
            self.field_data.append(datum)
            
class BedFile(object):
    def __init__(self):
        self.data = []

    def __iter__(self):
        for datum in self.data:
            yield datum

    def __len__(self):
        return len(self.data)

    def add_entry(self, entry):
        self.data.append(entry)
    def sort(self):
        self.data.sort(key = lambda a: a['start'])
    def cleanup(self):
        """
        Remove all duplicate entries, entries are mutable therefore not hashable
        on their own, but we can use the __repr__ of them to create a new list
        """
        unique = set()
        new_data = []
        for datum in self.data:
            uniq_str = datum.__repr__()
            if uniq_str in unique:
                continue
            else:
                unique.add(uniq_str)
                new_data.append(datum)
        self.data = new_data

    def from_bed_file(self, fname, header=0, names=None, dtypes=None):
        with open(fname, mode="r") as inf:
            for i in range(header):
                inf.readline()
            for line in inf:
                entry = BedEntry(names=names, dtypes=dtypes)
                entry.parse_from_line(line)
                self.data.append(entry)

    def write_bed_file(self, fname, header=None):
        with open(fname, mode="w") as outf:
            if header:
                for head_line in header:
                    outf.write(head_line)
            for entry in self:
                outf.write(str(entry))

    def to_gff(self, outfile):
        import gfftools
        with open(outfile, mode="w") as outf:
            for entry in self:
                gff_entry = gfftools.GffEntry()
                gff_entry.genome_name=entry["chrm"]
                gff_entry.data_origin = entry["name"]
                gff_entry.site_type = "."
                gff_entry.start = entry["start"] +1
                gff_entry.end = entry["end"]
                gff_entry.direction = entry["strand"]
                outf.write(str(gff_entry)+"\n")

if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    outfile = sys.argv[2]

    inbed = BedFile()
    inbed.from_bed_file(infile)
    inbed.to_gff(outfile)
