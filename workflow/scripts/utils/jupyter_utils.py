#!/usr/bin/env python3

class PDF(object):
  def __init__(self, pdf, size=(200,200)):
    self.pdf = pdf
    self.size = size

  def _repr_html_(self):
    return '<iframe src={0} width={1[0]} height={1[1]}></iframe>'.format(self.pdf, self.size)

  def _repr_latex_(self):
    return r'\includegraphics[width=1.0\textwidth]{{{0}}}'.format(self.pdf)
  
def get_ipysheet(cell_per_sample, sample):
  from ipysheet import sheet, column, to_dataframe
  import ipywidgets as w

  cell_list = cell_per_sample[sample]
  row_nb = len(cell_list)
  s = sheet(rows=row_nb, columns=2, column_headers=["cell", "selected?"])
  s.column_width=[8,2]
  s.layout = w.Layout(width='400px',height='100%')

  column(0, cell_list)
  column(1, [True]*row_nb)
  return s