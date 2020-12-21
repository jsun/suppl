# ruby tair_goslim2topgolist.rb > AT_GO_all.txt

require 'open-uri'
url = 'https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt'


atids = {}
goids = {}
atid2goid = {}
open(url) do |file|
  while line=file.gets
    x = line.chomp.split("\t")
    atid = x[0]
    goid = x[5]
    if atid =~ /^AT.G/ and goid =~ /^GO/
      atids[atid] = true
      goids[goid] = true
      atid2goid[atid] ||= []
      atid2goid[atid] << goid
    end
  end
end

atid2goid.sort.each do |atid, goids|
  puts [atid, goids.uniq.sort.join(",")].join("\t")
end



