_dict = {}
with open('Big5-ZhuYin.map', mode='r', encoding='big5hkscs') as fin:
    lines = fin.readlines()
    for line in lines:
        Big5, ZhuYins = line.split(' ', 1)
        ZhuYins = ZhuYins.replace('\n','').split('/')
        # 自己對到自己
        _dict[Big5] = [Big5]
        for ZhuYin in ZhuYins:
            if ZhuYin[0] in _dict.keys():
                _dict[ZhuYin[0]].append(Big5)
                _dict[ZhuYin[0]] = list(set(_dict[ZhuYin[0]]))
                _dict[ZhuYin[0]].sort()
            else:
                _dict[ZhuYin[0]] = [Big5]

with open( 'ZhuYin-Big5.map', mode='w', encoding='big5hkscs' ) as fout:
    for key in _dict.keys():
        fout.write(f'{key} {" ".join(_dict[key])}\n')
        # fout.write('{s1} {s2}\n'.format(s1=key, s2=" ".join(_dict[key])))