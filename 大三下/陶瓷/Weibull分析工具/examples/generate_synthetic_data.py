"""Generate explicitly synthetic strength data, not experiments or literature data."""
import csv
import random
from pathlib import Path

rng = random.Random(20261014)
output = Path(__file__).with_name('synthetic_strength.csv')
with output.open('w', newline='', encoding='utf-8') as stream:
    writer = csv.writer(stream)
    writer.writerow(['group', 'value', 'unit', 'data_kind'])
    for shape in (4, 12):
        for _ in range(60):
            writer.writerow([f'synthetic_m{shape}', f'{rng.weibullvariate(1200, shape):.6f}', 'MPa', 'synthetic'])
print(output)
