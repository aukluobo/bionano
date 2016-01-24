use strict;use warnings;

die "perl $0 <ctg.agp> <bionano.agp> <complete_relation.lst>" if @ARGV==0;

my ($ctg_agp,$bionano_agp,$cpt_rel)=@ARGV;

my %all;my $count=1;my %all_ctg;my %point;my %gap;
open IN,$bionano_agp or die;
my @pre=('CTG',0,0);
while(<IN>){
	if(/^#/){next;}
	chomp;
	my @aa=split;
	$all{$aa[0]}{$count}=[@aa];
	if(@aa>9){
		push @{$all_ctg{$aa[9]}{$aa[0]}},$count;print "$aa[0]\t$aa[9]\t$count\n";
		push @{$point{$aa[9]}{$aa[0]}},$aa[10];
		@pre=($aa[9],$aa[10],$pre[1]);
	}else{
		$gap{$pre[0]}{$pre[1]}=[$pre[2],$aa[5]];
	}
	$count+=1;
}
close IN;
$gap{$pre[0]}{$pre[1]}=[$pre[2],0];

my %ctg_agp;my %cpt_rel;my %all_bac;
read_ctg($ctg_agp,\%ctg_agp);

read_rel($cpt_rel,\%cpt_rel);
foreach my $ctg(keys %cpt_rel){
	foreach($cpt_rel{$ctg}){
		push @{$all_bac{$ctg}},$_;
	}
}

my %ex;my %link;
foreach my $sp_scf(keys %all){
	my $pre_ctg;my @sort;my @gap;my $pre_gap=0;my $gap_use=0;my $head_tail=0;
	my @sp_sort=sort {$a<=>$b} keys %{$all{$sp_scf}};
	#foreach my $sort(@sp_sort){
	for(my $i=0;$i<@sp_sort;$i++){
		my $sort=$sp_sort[$i];
		if(@{$all{$sp_scf}{$sort}}>9){
			unless($pre_ctg){
				$pre_ctg=$all{$sp_scf}{$sort}[9];
				push @sort,$all{$sp_scf}{$sort}[10];
				next;
			}
			if($pre_ctg ne $all{$sp_scf}{$sort}[9]){
				#deal with old CTG
				#foreach(@sort){$ex{$pre_ctg}{$_}=1;}
				#my $current_gap=$gap[-1];
				my $tmp_sort=1;my %tmp_sort;
				my %repeat;foreach(@sort){$repeat{$_}++;}
				foreach(keys %repeat){if($repeat{$_}>1){undef $repeat{$_};delete $repeat{$_};}}
				my $cha=$sort[0]-$sort[-1];
				for(my $i=0;$i<@sort;$i++){
					if($repeat{$sort[$i]}){
						my $check1=0;my $check2=0;
						if($i-1>0){$check1=$cha/($sort[$i-1]-$sort[$i]);}
						if($i+1<@sort){$check2=$cha/($sort[$i]-$sort[$i+1]);}
						if($check1 != 0){
							if($check2 != 0){
								if($check1/$check2<0){next;}
							}else{
								if($check1<0){next;}
							}
						}else{
							if($check2<0){next;}
						}
					}
					$tmp_sort{$sort[$i]}=$tmp_sort;
					$tmp_sort+=1;
				}
				@sort=sort {$tmp_sort{$a}<=>$tmp_sort{$b}} keys %tmp_sort;
				my $pre_c;my @pre_c;
				for(my $i=0;$i<@sort-1;$i++){
					unless($pre_c){$pre_c=$sort[$i+1]-$sort[$i];}
					unless($pre_c/($sort[$i+1]-$sort[$i])>0){push @pre_c,$i;}
					$pre_c=$sort[$i+1]-$sort[$i];
				}
				my $flag=0;
				if(@pre_c>0){
					print STDERR "$sp_scf bad BAC sort\n";
					$flag=1;
				}
				if($sort[0] eq (sort {$a<=>$b} keys %{$all{$sp_scf}})[0]){$head_tail=1;}
				foreach(@sort){print STDERR "$_\t";}print STDERR "\n";
				my @block=Find_Part(\%ctg_agp,\%cpt_rel,\%point,\%ex,\@sort,$pre_ctg,\%gap,$head_tail,$flag);
				push @{$link{$sp_scf}{$sp_sort[$i-2]}},$pre_ctg,@block;
				foreach my $ele(@block){
					print "$sp_scf\t$sp_sort[$i-2]\t$pre_ctg\t$ele\n";
				}

				$pre_ctg=$all{$sp_scf}{$sort}[9];
				$pre_gap=$gap[-1];
				@gap=();@sort=();
				push @sort,$all{$sp_scf}{$sort}[10];
			}else{
				$pre_ctg=$all{$sp_scf}{$sort}[9];push @sort,$all{$sp_scf}{$sort}[10];
			}
		}
		else{
			push @gap,$all{$sp_scf}{$sort}[5];
			$gap{$sort[-1]}=$all{$sp_scf}{$sort}[5];
		}
	}
	$head_tail=-1;
	my $tmp_sort=1;my %tmp_sort;
	my $cha=$sort[0]-$sort[-1];
	my %repeat;foreach(@sort){$repeat{$_}++;}foreach(keys %repeat){if($repeat{$_}>1){undef $repeat{$_};delete $repeat{$_};}}
	for(my $i=0;$i<@sort;$i++){
		if($repeat{$sort[$i]}){
			my $check1=0;my $check2=0;
			if($i-1>0){$check1=$cha/($sort[$i-1]-$sort[$i]);}
			if($i+1<@sort){$check2=$cha/($sort[$i]-$sort[$i+1]);}
			if($check1 != 0){
				if($check2 != 0){
					if($check1/$check2<0){next;}
				}else{
					if($check1<0){next;}
				}
			}else{
				if($check2<0){next;}
			}
		}
		$tmp_sort{$sort[$i]}=$tmp_sort;$tmp_sort+=1;
	}
	@sort=sort {$tmp_sort{$a}<=>$tmp_sort{$b}} keys %tmp_sort;
	my $pre_c;my @pre_c;
	for(my $i=0;$i<@sort-1;$i++){
		unless($pre_c){$pre_c=$sort[$i+1]-$sort[$i];}
		unless($pre_c/($sort[$i+1]-$sort[$i])>0){push @pre_c,$i;}
		$pre_c=$sort[$i+1]-$sort[$i];
	}
	my $flag=0;
	if(@pre_c>0){
		print STDERR "$sp_scf bad BAC sort\n";
		$flag=1;
	}
	my @block=Find_Part(\%ctg_agp,\%cpt_rel,\%point,\%ex,\@sort,$pre_ctg,\%gap,$head_tail,$flag);
	push @{$link{$sp_scf}{$sp_sort[-1]}},$pre_ctg,@block;
	foreach my $ele(@block){
		print "$sp_scf\t$sp_sort[-1]\t$pre_ctg\t$ele\n";
	}
}

#link the super_scaffold together
#first: find the head and tail linkage information

my %h_t;my %ht;
foreach my $sp_scf(keys %link){
	my @sp_sort =sort {$a<=>$b} keys %{$link{$sp_scf}};
	if(@sp_sort>1){
		my $head=$sp_sort[0];my $tail=$sp_sort[-1];
		my @head=@{$link{$sp_scf}{$head}};my @tail=@{$link{$sp_scf}{$tail}};
		my $head_ctg=shift @head;my $tail_ctg=shift @tail;
		print "$head_ctg\t$tail_ctg\n";
		$h_t{$head_ctg}{$sp_scf}{'head'}=[@head];
		$h_t{$tail_ctg}{$sp_scf}{'tail'}=[@tail];
		if(@head>1){
			$ht{$head_ctg}{$head[0]}{$sp_scf}='head1';
			$ht{$head_ctg}{$head[-1]}{$sp_scf}='head2';
		}else{
			$ht{$head_ctg}{$head[0]}{$sp_scf}='head0';
		}
		if(@tail>1){
			$ht{$tail_ctg}{$tail[0]}{$sp_scf}='tail1';
			$ht{$tail_ctg}{$tail[-1]}{$sp_scf}='tail2';
		}else{
			$ht{$tail_ctg}{$tail[0]}{$sp_scf}='tail0';
		}
		print "$head_ctg\t$sp_scf\thead aa\n$tail_ctg\t$sp_scf\ttail aa\n";
	}else{
		my @tmp=@{$link{$sp_scf}{$sp_sort[0]}};
		my $ctg=shift @tmp;
		$h_t{$ctg}{$sp_scf}{'both'}=[@tmp];
		if(@tmp>1){
			$ht{$ctg}{$tmp[0]}{$sp_scf}='both1';$ht{$ctg}{$tmp[-1]}{$sp_scf}='both2';
		}else{
			$ht{$ctg}{$tmp[0]}{$sp_scf}='both';
		}
		print "$ctg\t$sp_scf\tboth\n";
	}
}

my %anchor;my $anchor_sort=1;my %new_anchor;
foreach my $ctg(keys %h_t){
	my @sp_get=keys %{$h_t{$ctg}};
	if(@sp_get==1){
		#ctg on sp_scaffold
		print "$ctg\tonly on super_scaffold\n";
		$new_anchor{$ctg}{$sp_get[0]}=1;
		my $th=(keys %{$h_t{$ctg}{$sp_get[0]}})[0];
		#my $p1;
		if($th=~/head/){
			my $p1=$h_t{$ctg}{$sp_get[0]}{$th}[0];
			$anchor{$anchor_sort}{1}=[$sp_get[0],$ctg,$p1,$th];
			print "anchor_$anchor_sort\tnew_1\t$sp_get[0]\t$ctg\t$p1\t$th\n";$anchor_sort+=1;
		}
		elsif($th=~/tail/){
			my $p1=$h_t{$ctg}{$sp_get[0]}{$th}[-1];
			$anchor{$anchor_sort}{1}=[$sp_get[0],$ctg,$p1,$th];
			print "anchor_$anchor_sort\tnew_1\t$sp_get[0]\t$ctg\t$p1\t$th\n";$anchor_sort+=1;
		}else{
			my $p1=$h_t{$ctg}{$sp_get[0]}{$th}[0];
			my $p2=$h_t{$ctg}{$sp_get[0]}{$th}[-1];
			$anchor{$anchor_sort}{1}=[$sp_get[0],$ctg,$p1,$th];$anchor{$anchor_sort}{2}=[$sp_get[0],$ctg,$p2,$th];
			print "anchor_$anchor_sort\tnew_1\t$sp_get[0]\t$ctg\t$p1\t$th\n";
			print "anchor_$anchor_sort\tnew_1\t$sp_get[0]\t$ctg\t$p2\t$th\n";
			$anchor_sort+=1;
		}
		#$anchor{$anchor_sort}{1}=[$sp_get[0],$ctg,,$th];
		#print "anchor_$anchor_sort\tnew_1\t$_";
		#$anchor_sort+=1;
		print "\n";
	}elsif(@sp_get>=2){
		#multi-end
		my $cut=0;
		foreach my $spg(@sp_get){
			my @key=keys %{$h_t{$ctg}{$spg}};
			if(@key==1 and $key[0] ne 'both'){$cut+=1;}
		}
		if($cut>2){
			print "$ctg\t$cut\tcut CTG\n";
		}
		foreach my $bac_sort(sort {$a<=>$b} keys %{$ht{$ctg}}){
			my $sp_scf=(keys %{$ht{$ctg}{$bac_sort}})[0];
			print STDERR "$ctg\t$sp_scf\t$bac_sort\t$ht{$ctg}{$bac_sort}{$sp_scf}\n";
		}
		my @key=sort {$a<=>$b} keys %{$ht{$ctg}};
		my %new_sp;my $new_sp=1;my $last_link=0;my $last_add;
		for(my $i=0;$i<@key-1;$i++){
			my $p1=$key[$i];my $p2=$key[$i+1];
			my $p3;my $s3;my $marker_p3;my $p3_aft;my $p3_pre;
			if($i>0){
				$p3=$key[$i-1];
				$s3=(keys %{$ht{$ctg}{$p3}})[0];
				$marker_p3=$ht{$ctg}{$p3}{$s3};
				$p3_pre=substr($marker_p3,0,4);
				$p3_aft=$marker_p3;$p3_aft=~s/$p3_pre//;unless($p3_aft){$p3_aft=-1;}
			}else{$p3_aft=-1;}
			my $s1=(keys %{$ht{$ctg}{$p1}})[0];
			my $s2=(keys %{$ht{$ctg}{$p2}})[0];
			my $marker_p1=$ht{$ctg}{$p1}{$s1};
			my $marker_p2=$ht{$ctg}{$p2}{$s2};
			my $p1_pre=substr($marker_p1,0,4);my $p2_pre=substr($marker_p2,0,4);
			my $p1_aft=$marker_p1;$p1_aft=~s/$p1_pre//;
			my $p2_aft=$marker_p2;$p2_aft=~s/$p2_pre//;
			$last_add="$s2\t$ctg\t$p2\t$marker_p2\n";
			if($s1 eq $s2){
				push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
			}else{
				if($p1_pre eq $p2_pre){
					#my $p1_aft=$marker_p1;$p1_aft=~s/$p1_pre//;
					#my $p2_aft=$marker_p2;$p2_aft=~s/$p2_pre//;
					if($p1_pre eq 'head' and $p1_aft!=2 and $p2_aft!=2 or $p1_pre eq 'tail' and $p1_aft != 1 and $p2_aft !=1 or $p1_pre eq 'both'){
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
					}else{
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";
						$new_sp+=1;$last_link=0;
					}
				}else{
					if($p1_pre eq 'tail' and $p1_aft!=1){
						if($p2_pre eq 'both' or $p2_pre eq 'head' and $p2_aft !=2 and $p3_aft !=0 ){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
						}else{
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$new_sp+=1;$last_link=0;
						}
					}
					elsif($p1_pre eq 'head' and $p1_aft != 2){
						if($p2_pre eq 'both' or $p2_pre eq 'tail' and $p2_aft !=1 and $p3_aft !=0 ){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
						}else{
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$new_sp+=1;$last_link=0;
						}
					}elsif($p1_pre eq 'both'){
						if($p2_pre eq 'tail' and $p2_aft != 1 or $p2_pre eq 'head' and $p2_aft !=2){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
							print "$s1\t$ctg\t$p1\t$marker_p1\n";
						}else{
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$new_sp+=1;$last_link=0;
						}
					}
				}
			}
		}
		if($last_link){
			push @{$new_sp{$new_sp}},$last_add;
		}else{
			$new_sp+=1;push @{$new_sp{$new_sp}},$last_add;
		}
		foreach my $sot(sort {$a<=>$b} keys %new_sp){
			my $bb_sort=1;
			foreach(@{$new_sp{$sot}}){
				print "anchor_$anchor_sort\t$bb_sort\tnew_$sot\t$_";
				my @aa=split /\s+/,$_;
				$anchor{$anchor_sort}{$bb_sort}=[@aa];
				$bb_sort+=1;
			}
			$anchor_sort+=1;
			print "\n";
		}
	}
}

#merge sp_scf
my %mer;
foreach my $a_sort(sort {$a<=>$b} keys %anchor){
	my @b_sort=sort {$a<=>$b} keys %{$anchor{$a_sort}};
	if(@b_sort>1){
		my $st=$anchor{$a_sort}{$b_sort[0]}[0];
		my $ed=$anchor{$a_sort}{$b_sort[-1]}[0];
		$mer{$st}{$a_sort}='st';
		$mer{$ed}{$a_sort}='ed';
	}else{
		my $b_sort=$b_sort[0];
		my $type=$anchor{$a_sort}{$b_sort}[3];
		my $sp=$anchor{$a_sort}{$b_sort}[0];
		my $tp;
		if($type=~/head/){$tp='st';}elsif($type=~/tail/){$tp='ed'}
		$mer{$sp}{$a_sort}=$tp;
	}
}

my %sp_mer;
foreach my $sp(keys %mer){
	my @key=keys %{$mer{$sp}};
	if(@key==1){print "$sp\tsingle\n";}
	else{
		my ($srt1,$srt2)=@key;
		my @srt1=sort {$a<=>$b} keys %{$anchor{$srt1}};
		my @srt2=sort {$a<=>$b} keys %{$anchor{$srt2}};
		my @qu=(0,0);
		if($mer{$sp}{$srt1} eq 'ed'){$qu[0]=-1;}
		if($mer{$sp}{$srt2} eq 'ed'){$qu[0]=-1;}
		print "$sp\t$srt1\t$srt1[$qu[0]]\t$anchor{$srt1}{$srt1[$qu[0]]}[1]\tlink\t";
		print "$sp\t$srt2\t$srt2[$qu[0]]\t$anchor{$srt2}{$srt2[$qu[1]]}[1]\n";
		#push @{$sp_mer{$sp}},$anchor{$srt1}{$srt1[$qu[0]]}[1],$anchor{$srt2}{$srt2[$qu[1]]}[1];
		push @{$sp_mer{$sp}},$srt1,$srt1[$qu[0]],$srt2,$srt2[$qu[1]];
	}
}

#ouput from the start record  raw record on %link
my %use_sp;my %lone_ctg;my %final;my $final=1;
my %order;my %anchor_order;
foreach my $a_sort( sort {$a<=>$b} keys %anchor){
	my @b_sort=sort {$a<=>$b} keys %{$anchor{$a_sort}};
	my $use=0;
	foreach my $b_sort(@b_sort){
		my $sp=$anchor{$a_sort}{$b_sort}[0];
		#my $ctg=$anchor{$a_sort}{$b_sort}[1];

		if($use_sp{$sp}){$use=1;last;}
	}
	if($use){next;}
	my $head_sp=$anchor{$a_sort}{$b_sort[0]}[0];
	my $tail_sp=$anchor{$a_sort}{$b_sort[-1]}[0];
	my %tmp_use;
	my $tyh=$anchor{$a_sort}{$b_sort[0]}[3];
	my $tyt=$anchor{$a_sort}{$b_sort[-1]}[3];
	my $head_order=1;my $tail_order=1;
	if($tyh=~/head/){$head_order=-1;}
	if($tyt=~/tail/){$tail_order=-1;}
	my @head_path=Find_Path(1,$head_order,\%anchor,\%sp_mer,$a_sort,$b_sort[0],\%tmp_use);
	my @tail_path=Find_Path(-1,$tail_order,\%anchor,\%sp_mer,$a_sort,$b_sort[-1],\%tmp_use);
	#print "start end\n";
	my $len_head=@head_path;my $len_tail=@tail_path;
	my %last_sort;
	my $srt=1;
	print "head\t";
	for(my $i=0;$i<$len_head;$i++){
		$last_sort{$srt}=$head_path[$i];
		print "$srt\t$head_path[$i]\t";
		$srt+=1;
	}
	print "tail\t";
	for(my $i=$len_tail-2;$i>=0;$i--){
		$last_sort{$srt}=$tail_path[$i];
		print "$srt\t$tail_path[$i]\t";
		$srt+=1;
	}
	print "\n";
	#check sort 1 or -1
	my $plus=0;my $minus=0;
	my %tmp;
	my $st_sp;my $ed_sp;
	foreach my $srt(sort {$a<=>$b} keys %last_sort){
		my $asort=$last_sort{$srt};
		my @b_sort=sort {$a<=>$b} keys %{$anchor{$asort}};
		#my $start=$b_sort[0];my $end=$b_sort[-1];
		
		foreach my $b_sort(@b_sort){
			my $th=$anchor{$asort}{$b_sort}[3];
			my $spp=$anchor{$asort}{$b_sort}[0];
			unless($st_sp){$st_sp=$spp;}
			$ed_sp=$spp;
			push @{$tmp{$spp}},$th;
		}
	}
	#my %order;
	foreach my $sppp(keys %tmp){
		if(@{$tmp{$sppp}}<2){
			if($tmp{$sppp}[0]=~/both/){
				next;
			}else{
				if($tmp{$sppp}[0]=~/head0/){
					if($sppp eq $st_sp){$minus+=1;$order{$sppp}=-1;}
					if($sppp eq $ed_sp){$plus+=1;$order{$sppp}=1;}
				}elsif($tmp{$sppp}[0]=~/tail0/){
					if($sppp eq $st_sp){$plus+=1;;$order{$sppp}=1;}
					if($sppp eq $ed_sp){$minus+=1;$order{$sppp}=-1;}
				}
			}
		}else{
			my $aa=$tmp{$sppp}[0];my $bb=$tmp{$sppp}[1];
			if($aa=~/head/){
				if($aa=~/2$/){$minus+=1;;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;;$order{$sppp}=1;}
			}elsif($aa=~/tail/){
				if($aa=~/2$/){$minus+=1;;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;;$order{$sppp}=1;}
			}elsif($aa=~/both/){
				if($aa=~/2$/){$minus+=1;;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;;$order{$sppp}=1;}
			}
		}
	}
	my $order;
	if($minus>$plus){$order=-1;}
	else{
		$order=1;
	}
	my @srt=sort {$a<=>$b} keys %last_sort;
	if($order<0){@srt=sort {$b<=>$a} @srt;}
	#push @{$final{$final}},@srt,$order;$final+=1;
	for(my $i=0;$i<@srt;$i++){
		my $asort=$last_sort{$srt[$i]};
		push @{$final{$final}},$asort;
		foreach my $bb_sort(keys %{$anchor{$asort}}){
			my $ctg=$anchor{$asort}{$bb_sort}[1];
			$lone_ctg{$ctg}{$asort}=$bb_sort;
			my $sp=$anchor{$asort}{$bb_sort}[0];
			$use_sp{$sp}=1;
		}
	}
	push @{$final{$final}},$order;
	$final+=1;
}

my %no_split;
foreach my $ctg(keys %lone_ctg){
	my @srt=keys %{$lone_ctg{$ctg}};
	if(@srt==1){
		$no_split{$ctg}=1;
		print "$ctg\tno split\n";
	}
}
%ex=();
my %bac_sort;my %on_bac;
foreach my $ff(sort {(keys $final{$b})<=>(keys $final{$a})} keys %final){
	my $bac_sort=1;
	my @asort=@{$final{$ff}};
	my $ff_order=pop @asort;
	my $out=join "\t",@asort;
	#print "final sort\t$ff\t$order\t$out\n";
	my $po=1;
	foreach my $ii(@asort){
		my ($asort,$order)=split /_/,$ii;
		if($ex{$asort}){$po=0;last;}
		$ex{$asort}=1;
	}
	print "final sort\t$ff\t$ff_order\t$out\n" if $po;
	if($po){
		#this sort will be used
		my $plus=0;my $minus=0;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			if($asort_order>0){$plus+=1;}if($asort_order<0){$minus+=1;}
		}
		my $multi=1;
		if($minus > $plus){@asort=reverse @asort;$multi=-1;}
		my %last_sp;my $cu=1;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			$asort_order*=$multi;
			my @bsort=sort {$a<=>$b} keys %{$anchor{$asort}};
			if($asort_order<0){@bsort=reverse @bsort;}
			foreach my $bsort(@bsort){
				my @tmp=@{$anchor{$asort}{$bsort}};
				$last_sp{$tmp[0]}{'sort'}=$cu;
				$last_sp{$tmp[0]}{$asort}=$asort_order;
				$cu+=1;
			}
		}
		foreach my $sp(sort {$last_sp{$a}{'sort'}<=>$last_sp{$b}{'sort'}} keys %last_sp){
			my $sp_order;
			#if($order{$sp}){$sp_order=$order{$sp};}else{$sp_order=1;}
			#my @abc=@{$last_sp{$sp}};print "finalCTG\t$ff\t$sp\t",join "\t",@abc,"\n";
			my $asort;
			foreach my $aas(keys %{$last_sp{$sp}}){
				if($aas=~/sort/){next;}
				$asort.="|$aas";
				unless($sp_order){$sp_order=$last_sp{$sp}{$aas};}else{$sp_order*=$last_sp{$sp}{$aas};}
			}
			print "finalCTG\t$ff\t$sp\t$sp_order\n";
			#all sp informaiton are stored in %link
			my @key=keys %{$link{$sp}};
			if($sp_order>0){@key=sort {$a<=>$b} @key;}
			if($sp_order<0){@key=sort {$b<=>$a} @key;}
			foreach my $key(@key){
				my @key1=@{$link{$sp}{$key}};
				my $ctg=shift @key1;
				#if($sp_order>0){@key1=sort {$a<=>$b} @key1;}
				if($sp_order<0){@key1=sort {$b<=>$a} @key1;}
				foreach(@key1){
					print "finalCTG\t$ff\t$asort\t$sp\t$ctg\t$_\n";
					push @{$bac_sort{$ff}{$bac_sort}},[$ctg,$_];
					$on_bac{$ctg}{$_}=[$sp,];
					$bac_sort+=1;
				}
			}
		}
	}
}

#find the missing BAC sort
my $input
foreach my $ctg(keys %cpt_rel){
	my %miss;my %out_on;
	foreach my $ss(sort {$a<=>$b} keys %{$cpt_rel{$ctg}}){
		push @{$out_on{$ss}},$ctg;
		if($on_bac{$ctg} and $on_bac{$ctg}{$ss}){
			push @{$out_on{$ss}},"on";
			next;
		}
		push @{$out_on{$ss}},"out";
		$miss{$ss}=1;
		print "$ctg\t$ss\tmiss\n";
	}
	#find a place to locate the miss bac sort and output the last bac sort;some place may need to split by the overlap
	#information because they can not link together
	if($on_bac{$ctg} and keys %miss >0){
		my @key=keys %out_on;my $st=-1;my $ed=0;
		for(my $i=0;$i<keys @key;$i+=1){
			if($out_on{$key[$i]}[1] eq "out"){
				if($st<0){$st=$i;}
			}
			if($out_on{$key[$i]}[1] eq "on" and $st>=0){
				unless ($ed) {
					$ed=$i;
				}
			}
			if($st>=0 and $ed>=0){
				#get out block; start to decide which block it belong to
				my $pre_on=$st-1;if($st==0){$pre_on=-1;}
				my $aft_on=$ed+1;if($ed+1>@key-1){$aft_on=-1;}
				if($pre_on<0){
					#no pre
					if($aft_on<0){
						#no anchor
						print "input\t$out_on{$key[$i]}[0]\tun_anchor\n";
					}else{
						#input on aft
					}
				}else{
					if($aft_on<0){
						#input on pre
					}else{
						#both. check which one is properble.
					}
				}
			}
		}
	}
	
}











#sub routine

sub Find_Path
{
	my ($type,$order,$hash,$mer_hash,$asort,$bsort,$use)=@_;
	my $pre_sp=$hash->{$asort}{$bsort}[0];
	#print "start $pre_sp\t$asort\t$bsort\n";
	my $nextasort;my $nextbsort;
	my @back;
	$use->{$asort}=1;
	my $raw_order=$order;
	if($mer_hash->{$pre_sp}){
		if($asort == $mer_hash->{$pre_sp}[0] and $bsort == $mer_hash->{$pre_sp}[1]){
			$nextasort=$mer_hash->{$pre_sp}[2];$nextbsort=$mer_hash->{$pre_sp}[3];
		}else{
			$nextasort=$mer_hash->{$pre_sp}[0];$nextbsort=$mer_hash->{$pre_sp}[1];
		}
		my @next=sort {$a<=>$b} keys %{$hash->{$nextasort}};
		my $nextb;
		if($nextbsort == $next[0]){
			$nextb=$next[-1];
			if($type==1){$order=-1;}
			elsif($type==-1){$order=1;}
		}else{
			$nextb=$next[0];
			if($type==1){$order=1;}
			elsif($type==-1){$order=-1;}
		}
		if(@next==1){$order=1;$use->{$nextasort}+=1;}
		if($use->{$nextasort}){
			;
		}else{
			#$tht=$hash->{$nextasort}{$nextb}[3];
			@back=Find_Path($type,$order,$hash,$mer_hash,$nextasort,$nextb,$use);
		}
	}
	my $out="$asort";if($raw_order){$out.="_$raw_order";}
	push @back,$out;

	return @back;
}


sub Find_Part
{
	my ($ctgagp,$cptrel,$point,$ex,$sort,$ctg,$gap,$head_tail,$flag)=@_;
	print STDERR "checking $ctg\n";
	my @ss=@{$sort};
	my @block;
	if($flag){
		#find the bad sort point
		my %keep;
		for(my $i=0;$i<@ss-2;$i+=1){
			my $dvp=$ss[$i+1]-$ss[$i];my $dva=$ss[$i+2]-$ss[$i+1];
			if($dvp/$dva<0){
				$keep{$i}=1;$keep{$i+1}=1;
			}
		}
		my %mul;my $tt=1;
		for(my $i=0;$i<@ss;$i+=1){
			if($keep{$i}){$tt+=1;next;}
			push @{$mul{$tt}},$ss[$i];
		}
		my @tt=sort {$a<=>$b} keys %mul;
		if(@{$mul{$tt[-1]}}==2){pop @tt;}
		foreach my $t(@tt){
			if(@{$mul{$t}}<=1){next;}
			push @block,[@{$mul{$t}}];
		}
	}else{
		if(@ss>1){
			push @block,[@ss];
		}
	}
	
	#check sort is bloken or not
	my @output;
	foreach my $ss(@block){
		#if(@{$ss}==1){push @output,@{$ss};next;}
		@ss=sort {$a<=>$b} @{$ss};
		my %bloken;
		
		foreach my $sp(keys %{$point->{$ctg}}){
			my @tmp_ss=sort {$a<=>$b} @{$point->{$ctg}{$sp}};
			if($tmp_ss[0]==$ss[0]){next;}
			my %tmp;
			foreach my $rss(@ss){$tmp{$rss}='r';}
			foreach my $tss(@tmp_ss){$tmp{$tss}='t';}
			#my $chr;
			my %cc;
			my @cc=sort {$a<=>$b} keys %tmp;
			for(my $i=0;$i<@cc-1;$i++){
				if($tmp{$cc[$i]} ne $tmp{$cc[$i+1]}){
					$cc{"$cc[$i]\_$cc[$i+1]"}=1;
				}
			}
			my @key=keys %cc;
			if(@key>1){
				#bloken
				foreach my $cc_key(keys %cc){
					my @cc_key=split /_/,$cc_key;
					if($tmp{$cc_key[0]} eq 'r'){push @{$bloken{$tmp{$cc_key[0]}}},$cc_key[1];}
					if($tmp{$cc_key[1]} eq 'r'){push @{$bloken{$tmp{$cc_key[1]}}},$cc_key[0];}
				}
			}
		}
		#my @output;
		if(keys %bloken >0){
			for(my $i=0;$i<@{$sort}-1;$i++){
				my $pre=$sort->[$i];my $aft=$sort->[$i];

				if($bloken{$pre} and $bloken{$aft}){
					my $gap=$gap->{$ctg}{$pre}[1];
					my @get=(0,-1);if($pre > $aft){@get=(-1,0);}
					my $get_pre=(sort {$a<=>$b} @{$bloken{$pre}})[$get[0]];
					my $get_aft=(sort {$a<=>$b} @{$bloken{$aft}})[$get[1]];
					print STDERR "double check $pre\t$get_pre and $aft\t$get_aft\n";
					my @check_pre=Check_Bpoint($pre,$get_pre,$ctgagp,$cptrel);
					my @check_aft=Check_Bpoint($aft,$get_aft,$ctgagp,$cptrel);
					push @output,$pre,@check_pre,@check_aft,$aft;
				}else{
					my @tmp_output=($pre..$aft);
					if($pre>$aft){@tmp_output=($aft..$pre);@tmp_output=reverse @tmp_output;}
					push @output,@tmp_output;
				}
			}
		}else{
			my $aa=$sort->[0];my $bb=$sort->[-1];
			my @tmp_output=($aa..$bb);
			if($aa>$bb){@tmp_output=($bb..$aa);@tmp_output=reverse @tmp_output;}
			push @output,@tmp_output;
		}
	}
	#if($flag){@output=@{$sort};}
	my %output;
	if(@output){
		print STDERR "test\t",join "\t",@output,"\n";
		my ($st_o,$ed_o)=(sort {$a<=>$b} @output)[0,-1];
		my $ou=0;my @tmp=@output;@output=();
		foreach my $ele(@{$sort}){
			if($ele == $st_o or $ele== $ed_o){
				$ou+=1;
				if($ou==1){
					foreach(@tmp){
						push @output,$_;
					}
					$ou+=1;
				}
			}elsif($ou == 0 or $ou == 3){
				push @output,$ele;
			}
		}
	}else{
		push @output,@{$sort};
	}
	foreach(@output){print STDERR "$_\t";}
	print STDERR "\n";
	foreach(@{$sort}){print STDERR "$_\t";}
	print STDERR "\n";
	return @output;
	#extend head or tail
	if($head_tail==1){
		#extend head block
		#my $add=1;if($output[0]-$output[1]>0){$add=-1;}
		my $st=$output[0];
		while($st>1){
			$st-=1;
			#
		}
	}
}

sub Check_Bpoint
{
	my ($pre,$gpre,$ctgagp,$cptrel)=@_;
	if(abs($pre-$gpre)==1){
		return $pre;
	}else{
		#
	}
}

sub read_ctg
{
	my ($file,$hash)=@_;
	open FI,$file or die $!;
	my $si=1;
	while(<FI>){
		if(/^#/){next;}
		chomp;
		my @aa=split;
		$hash->{$aa[0]}{$si}=[@aa];
		$si++;
	}
	close FI;
}

sub read_rel
{
	my ($file,$hash)=@_;
	open FI,$file or die $!;
	while(<FI>){
		if(/bac/){
			my @aa=split;
			if($aa[-1]=~/_(bac.*)_(\d+)$/){
				my $ctg=$`;my $bac=$1;my $ss=$2;
				my $length=0;
				if(@aa>4){$length=$aa[1]-$aa[3];}else{$length=$aa[1];}
				$hash->{$ctg}{$ss}=[$bac,$length];
			}
		}
	}
	close FI;
}


