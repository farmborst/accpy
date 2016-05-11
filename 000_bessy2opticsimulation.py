# -*- coding: utf-8 -*-
''' accpy example 0
author:     felix.kramer(at)physik.hu-berlin.de
'''

from accpy.simulate.lsd import lsd


bessy2lattices = ['bessy2injectionline',
                  'bessy2booster',
                  'bessy2transfer',
                  'bessy2ring']
lsd(bessy2lattices[2], slic=int(1e3), save=False, ft='pdf',
    plotstandard='presentation_1920x1080', scale=[1, 1])
