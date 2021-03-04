import os
import sys



N = 450000000




def get_intermediate_position(gff_fpath):
    
    inter_pos = {}
    
    with open(gff_fpath, 'r') as infh:
        for buf in infh:
            if buf[0] == '#':
                continue
            
            chr_name, chr_ver, feature, pos_start, pos_end, \
                    qual, strand, status, attr = buf.replace('\n', '').split('\t')
            
            pos_start = int(pos_start)
            pos_end = int(pos_end)
            
            if chr_name not in inter_pos:
                inter_pos[chr_name] = {'inter_left': None, 'inter_right': None, 'inter': None}
            
            if pos_start < N:
                if inter_pos[chr_name]['inter_left'] is None or \
                    inter_pos[chr_name]['inter_left'] < pos_start:
                    inter_pos[chr_name]['inter_left'] = pos_start
                
            if pos_end > N:
                if inter_pos[chr_name]['inter_right'] is None or \
                    pos_end < inter_pos[chr_name]['inter_right']:
                    inter_pos[chr_name]['inter_right'] = pos_end
    
    print(inter_pos)
    
    # claculate inter mediate position
    inter = {}
    for k in inter_pos.keys():
        if inter_pos[k]['inter_right'] is not None:
            inter[k] = int((inter_pos[k]['inter_left'] + inter_pos[k]['inter_right']) / 2)
        else:
            inter[k] = None
            
    
    return inter





def write_gff(inter_pos, in_gff_fpath, out_gff_fpath):
    
    with open(in_gff_fpath, 'r') as infh, open(out_gff_fpath, 'w') as outfh:
        for buf in infh:
            if buf[0] == '#':
                continue
            
            chr_name, chr_ver, feature, pos_start, pos_end, \
                    qual, strand, status, attr = buf.replace('\n', '').split('\t')
            pos_start = int(pos_start)
            pos_end = int(pos_end)
            
            if inter_pos[chr_name] is None:
                part = 'part1'
            elif pos_start < inter_pos[chr_name] and pos_end < inter_pos[chr_name]:
                part = 'part1'
            elif inter_pos[chr_name] < pos_start and inter_pos[chr_name] < pos_end:
                part = 'part2'
                pos_start = pos_start - inter_pos[chr_name]
                pos_end = pos_end - inter_pos[chr_name]
            else:
                raise ValueError()
            
            new_chr_name = chr_name[0:5] + '_' + part
            new_gff_records = '\t'.join([new_chr_name, chr_ver, feature, str(pos_start), str(pos_end),
                                            qual, strand, status, attr])
            outfh.write(new_gff_records + '\n')





def __read_seq(fa_fpath):
    
    entry_id = ''
    entry_seq = ''
    
    with open(fa_fpath, 'r') as infh:
        for buf in infh:
            buf = buf.replace('\n', '')
            if buf[0] == '>':
                if entry_id != '' and entry_seq != '':
                    yield entry_id, entry_seq
                
                entry_id = buf[1:]
                entry_seq = ''
            else:
                entry_seq = entry_seq + buf
    
    yield entry_id, entry_seq




def __write_seq(outfh, seq):
    reformed_seq = ''
    
    i = 0
    n = 100
    while i <= len(seq):
        if i % 1e8 == 1e8 - 1:
            print(i)
        outfh.write(seq[i:(i+n)] + '\n')
        i = i + n
    




def write_seq(outfh, entry_id, entry_seq, inter_pos):
    print(entry_id)
    
    if inter_pos[entry_id] is not None:
        # part 1
        entry_id_part1 = entry_id[0:5] + '_part1'
        entry_seq_part1 = entry_seq[0:(inter_pos[entry_id])]
        outfh.write('>{}\n'.format(entry_id_part1))
        __write_seq(outfh, entry_seq_part1)
    
        # part 2
        entry_id_part2 = entry_id[0:5] + '_part2'
        entry_seq_part2 = entry_seq[(inter_pos[entry_id]):]
        outfh.write('>{}\n'.format(entry_id_part2))
        __write_seq(outfh, entry_seq_part2)
    
    else:
        entry_id_part1 = entry_id[0:5] + '_part1'
        entry_seq_part1 = entry_seq
        outfh.write('>{}\n'.format(entry_id_part1))
        __write_seq(outfh, entry_seq_part1)
        


    
def write_fasta(inter_pos, in_fa_fpath, out_fa_fpath):
    with open(out_fa_fpath, 'w') as outfh:
        for entry_id, entry_seq in __read_seq(in_fa_fpath):
            write_seq(outfh, entry_id, entry_seq, inter_pos)

    



if __name__ == '__main__':
    
    if len(sys.argv) != 5:
        raise ValueError('Err!')
    
    input_fasta_fpath = sys.argv[1]
    input_gff_fpath = sys.argv[2]
    output_fasta_fpath = sys.argv[3]
    output_gff_fpath = sys.argv[4]
    
    
    inter_pos = get_intermediate_position(input_gff_fpath)
    write_gff(inter_pos, input_gff_fpath, output_gff_fpath)
    write_fasta(inter_pos, input_fasta_fpath, output_fasta_fpath)
    

