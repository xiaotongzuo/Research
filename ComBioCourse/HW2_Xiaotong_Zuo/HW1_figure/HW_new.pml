bg white
load Downloads/1yy8.pdb
select L, resi 1-107
show cartoon, L
orient L
select L, resi 1-107 and chain A
orient L
hide line
cmd.spectrum("count",selection="(L)&elem C")
cmd.disable('L')
cmd.enable('L')
cmd.disable('1yy8')
cmd.enable('1yy8',1)
cmd.disable('L')
cmd.enable('L')
select L, resi 1-107 and chain A
select L, resi 1-107 + chain A
select L, resi 1-107 and chain A
cmd.select('L',"byresi((((L) or byresi((1yy8`2156))) and not ((byresi((1yy8`2156))) and byresi(L))))",enable=1)
cmd.select('L',"byresi((((L) or byresi((1yy8`2156))) and not ((byresi((1yy8`2156))) and byresi(L))))",enable=1)
hide all
show cartoon, L
cmd.disable('L')
select F21, chain A and resi 21
select L73, resi 73 and chain A
show stick, F21 + L73
wizard distance
refresh_wizard
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.select('sele','none')
cmd.select('sele',"((((sele) or ((1yy8`543))) and not ((((1yy8`543))) and (sele))))",enable=1)
cmd.get_wizard().do_select('''sele''')
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.select('_indicate_mw',"((((_indicate_mw) or ((1yy8`145))) and not ((((1yy8`145))) and (_indicate_mw))))",enable=1)
cmd.get_wizard().do_select('''_indicate_mw''')
cmd.get_wizard().do_dirty()
cmd.get_wizard().do_dirty()
cmd.set_wizard()
space cmyk
ray 1200, 800
cmd.log_close()
