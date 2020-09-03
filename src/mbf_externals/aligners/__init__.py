from .star import STAR
from .subread import Subread, Subjunc
from .bowtie import Bowtie
from .salmon import Salmon
from .bwa import BWA

all = [Bowtie, Subread, STAR, Salmon, Subjunc, BWA]
