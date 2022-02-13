BEGIN{
  Path=Folder"/Chrom_Sep"
}
{
  print $0 > Path"/"$1".sep"
}
