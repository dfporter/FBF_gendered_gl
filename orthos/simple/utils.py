import re


def label_cleaner(x):
    x = re.sub('Fbf1', 'FBF-1', x)
    x = re.sub('Fbf2', 'FBF-2', x)
    x = re.sub('Fbf', 'FBF', x)
    x = re.sub('Block', 'Block ', x)
    x = re.sub('_', ' ', x)
    x = re.sub('20C', '20°', x)
    x = re.sub('25C', '25°', x)
    x = re.sub('OO', ' OO ', x)
    x = re.sub('SP', ' SP ', x)
    x = re.sub('Top500', 'Top 500 ', x)
    return x