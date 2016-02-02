use strict;use warnings;

die "perl $0 <ctg.agp> <bionano.agp> <complete_relation.lst>" if @ARGV==0;

my ($ctg_agp,$bionano_agp,$cpt_rel)=@ARGV;

my %allp;my $count=1;my %all_ctg;my %point;my %gap;
open IN,$bionano_agp or die;
my @pre=('CTG',0,0);
my %bac_num;my %same_bac;
while(<IN>){
	if(/^#/){next;}
	chomp;
	my @aa=split;
	$allp{$aa[0]}{$count}=[@aa];
	if(@aa>9){
		push @{$all_ctg{$aa[9]}{$aa[0]}},$count;print "$aa[0]\t$aa[9]\t$count\n";
		push @{$point{$aa[9]}{$aa[0]}},$aa[10];
		@pre=($aa[9],$aa[10],$pre[1]);
		my $bac=$aa[5];$bac=~s/\d+$//;
		$bac_num{$aa[0]}{$bac}=1;
		push @{$same_bac{$aa[0]}{$bac}},$count;
	}else{
		$gap{$pre[0]}{$pre[1]}=[$pre[2],$aa[5]];
	}
	$count+=1;
}
close IN;

my %lone_ctg_sp;
foreach my $ctg(keys %all_ctg){
	my @sp=keys %{$all_ctg{$ctg}};
	if(@sp==1){
		$lone_ctg_sp{$ctg}=1;
	}
}

$gap{$pre[0]}{$pre[1]}=[$pre[2],0];

my %ctg_agp;my %cpt_rel;my %all_bac;
read_ctg($ctg_agp,\%ctg_agp);

read_rel($cpt_rel,\%cpt_rel);
foreach my $ctg(keys %cpt_rel){
	foreach($cpt_rel{$ctg}){
		push @{$all_bac{$ctg}},$_;
	}
}

#discard one sp contain one BAC
foreach my $sp(keys %bac_num){
	my @key=keys %{$bac_num{$sp}};
	if(@key==1){
		print STDERR "discard\t$sp\t$key[0]\n";
		undef $allp{$sp};delete $allp{$sp};
		undef $same_bac{$sp};delete $same_bac{$sp};
	}
}

my %ban;
foreach my $sp(keys %same_bac){
	my @cu=sort {$a<=>$b} keys %{$allp{$sp}};

	foreach my $bac(keys %{$same_bac{$sp}}){
		my @tmp=@{$same_bac{$sp}{$bac}};
		if(@tmp==1){next;}
		print STDERR "$bac\t",join "\t",@tmp,"\n";
		for(my $i=1;$i<@tmp;$i+=1){
			if(abs($tmp[$i-1]-$tmp[$i])==2){
				#discard pre one
				print STDERR "same\t$sp\t$bac\t$tmp[$i-1]\t$tmp[$i]\n";
				if($tmp[$i-1] == $cu[0]){
					$ban{$tmp[$i-1]}=1;my $nx=$tmp[$i-1]+1;$ban{$nx}=1;
					print STDERR "ban\t$sp\t$bac\t$tmp[$i-1]\t$nx\n";
				}
				if($tmp[$i]==$cu[-1]){
					$ban{$tmp[$i]}=1;my $nx=$tmp[$i]-1;$ban{$nx}=1;
					print STDERR "ban\t$sp\t$bac\t$tmp[$i]\t$nx\n";
				}
			}
		}
	}
}
my %all;
foreach my $sp(keys %allp){
	my @cu=sort {$a<=>$b} keys %{$allp{$sp}};
	foreach my $cu(@cu){
		if($ban{$cu}){next;}
		$all{$sp}{$cu}=$allp{$sp}{$cu};
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
					print "link\t$sp_scf\t$sp_sort[$i-2]\t$pre_ctg\t$ele\n";
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
		print "link\t$sp_scf\t$sp_sort[-1]\t$pre_ctg\t$ele\n";
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
		my %new_sp;my $new_sp=1;my $last_link=0;my $last_add;my %on;
		
		#make sure that the both is neibor
		my %tmp;
		for(my $i=0;$i<@key;$i++){
			my $p1=$key[$i];
			my $s1=(keys %{$ht{$ctg}{$p1}})[0];
			my $marker_p1=$ht{$ctg}{$p1}{$s1};
			if($marker_p1=~/both/){
				push @{$tmp{$s1}{$marker_p1}},$i;
			}
		}
		my %discard;my %to;
		foreach my $sp(keys %tmp){
			my @key1=keys %{$tmp{$sp}};
			if(@key1 == 2 and $key1[0]=~/both/){
				my @key2=sort {$a<=>$b} ($tmp{$sp}{$key1[0]}[0],$tmp{$sp}{$key1[1]}[0]);
				if(abs($key2[0]-$key2[1])==1){
					next;
				}else{
					if(abs($key2[0]-$key2[1])>2){
						#move them together
						$to{$key2[0]}=$key2[1];
						$discard{$key2[1]}=1;
					}else{
						#discard the middle one
						my $p1=$key[$key2[1]-1];
						my $s1=(keys %{$ht{$ctg}{$p1}})[0];
						my $marker_p1=$ht{$ctg}{$p1}{$s1};
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=0;
						print "outfirst\t$s1\t$ctg\t$p1\t$marker_p1\n";
						$new_sp+=1;
						$discard{$key2[1]-1}=1;
					}
				}
			}
		}
		my @key_tmp;my %bk;my $j=0;
		for(my $i=0;$i<@key;$i++){
			if($discard{$i}){print "discard\t$i\t$key[$i]\n";next;}
			if($to{$i}){
				$key_tmp[$j]=$key[$i];
				$key_tmp[$j+1]=$key[$to{$i}];$bk{$j+1}=1;$j+=2;
				print "fow\t$i\t$key[$i]\n";
				next;
			}
			$key_tmp[$j]=$key[$i];$j+=1;
		}

		my @key_old=@key;
		@key=@key_tmp;
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
			if($bk{$i}){
				push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=0;
				print "distest\t$s1\t$ctg\t$p1\t$marker_p1\tdis\n";
				$new_sp+=1;next;
			}
			if($s1 eq $s2){
				push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
			}else{
				if($p1_pre eq $p2_pre){
					#my $p1_aft=$marker_p1;$p1_aft=~s/$p1_pre//;
					#my $p2_aft=$marker_p2;$p2_aft=~s/$p2_pre//;
					my $abc=$i-1;
					if($on{$abc} and $on{$abc}==1 and $p1_pre ne 'both'){
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";
						print "$s1\t$ctg\t$p1\t$marker_p1\n";
						$new_sp+=1;$last_link=0;
						next;
					}
					if($p1_pre eq 'head' and $p1_aft!=2 and $p2_aft!=2 or $p1_pre eq 'tail' and $p1_aft != 1 and $p2_aft !=1 or $p1_pre eq 'both'){
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;
						$on{$i}=1;
						print "$s1\t$ctg\t$p1\t$marker_p1\tkk\n";
					}else{
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";
						$new_sp+=1;$last_link=0;
					}
				}else{
					my $abc=$i-1;
					if($on{$abc} and $on{$abc}==1 and $p1_aft == 0){
						push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";
						print "$s1\t$ctg\t$p1\t$marker_p1\n";
						$new_sp+=1;$last_link=0;
						next;
					}
					if($p1_pre eq 'tail' and $p1_aft!=1){
						if($p2_pre eq 'both' or $p2_pre eq 'head' and $p2_aft !=2 and $p3_aft !=0 ){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;$on{$i}=1;
						}else{
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$new_sp+=1;$last_link=0;
						}
					}
					elsif($p1_pre eq 'head' and $p1_aft != 2){
						if($p2_pre eq 'both' or $p2_pre eq 'tail' and $p2_aft !=1 and $p3_aft !=0 ){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;$on{$i}=1;
						}else{
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$new_sp+=1;$last_link=0;
						}
					}elsif($p1_pre eq 'both'){
						if($p2_pre eq 'tail' and $p2_aft != 1 or $p2_pre eq 'head' and $p2_aft !=2){
							push @{$new_sp{$new_sp}},"$s1\t$ctg\t$p1\t$marker_p1\n";$last_link=1;$on{$i}=1;
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
			print "anchor_$anchor_sort\t$ctg\txxx\t$sot\n";
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
		if($mer{$sp}{$srt2} eq 'ed'){$qu[1]=-1;}
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
	if($tyh=~/head/){$head_order=1;}
	if($tyt=~/tail/){$tail_order=1;}
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
		my $asort=(split /_/,$last_sort{$srt})[0];
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
		print "check $sppp\n";
		if(@{$tmp{$sppp}}<2){
			if($tmp{$sppp}[0]=~/both/){
				next;
			}else{
				if($tmp{$sppp}[0]=~/head0/){
					if($sppp eq $st_sp){$minus+=1;$order{$sppp}=-1;}
					if($sppp eq $ed_sp){$plus+=1;$order{$sppp}=1;}
				}elsif($tmp{$sppp}[0]=~/tail0/){
					if($sppp eq $st_sp){$plus+=1;$order{$sppp}=1;}
					if($sppp eq $ed_sp){$minus+=1;$order{$sppp}=-1;}
				}
			}
		}else{
			my $aa=$tmp{$sppp}[0];my $bb=$tmp{$sppp}[1];
			if($aa=~/head/){
				if($aa=~/2$/){$minus+=1;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;$order{$sppp}=1;}
			}elsif($aa=~/tail/){
				if($aa=~/2$/){$minus+=1;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;$order{$sppp}=1;}
			}elsif($aa=~/both/){
				if($aa=~/2$/){$minus+=1;$order{$sppp}=-1;}
				if($aa=~/1$/){$plus+=1;$order{$sppp}=1;}
			}
		}
	}
	my $order;
	if($minus>$plus){$order=-1;}
	else{
		$order=1;
	}
	my @srt=sort {$a<=>$b} keys %last_sort;
	#if($order<0){@srt=sort {$b<=>$a} @srt;}
	#push @{$final{$final}},@srt,$order;$final+=1;
	print "hehe\t";
	for(my $i=0;$i<@srt;$i++){
		my $asort=$last_sort{$srt[$i]};
		print "$asort\t";
		push @{$final{$final}},$asort;
		my ($aa_sort,$tmp_order)=split /_/,$asort;
		foreach my $bb_sort(keys %{$anchor{$aa_sort}}){
			my $ctg=$anchor{$aa_sort}{$bb_sort}[1];
			$lone_ctg{$ctg}{$aa_sort}=$bb_sort;
			my $sp=$anchor{$aa_sort}{$bb_sort}[0];
			$use_sp{$sp}=1;
			print "$sp use\t";
		}
	}
	print "\n";
	$order=1;
	push @{$final{$final}},$order;
	$final+=1;
}

#foreach my $ssp(keys %order){print "checkorder\t$ssp\t$order{$ssp}\n";}

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
foreach my $ff(sort {@{$final{$b}}<=>@{$final{$a}}} keys %final){
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
	print "final sort\t$ff\t$ff_order\t$out\ta\n" if $po;
	if($po){
		#this sort will be used
		my $plus=0;my $minus=0;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			if($asort_order>0){$plus+=1;}if($asort_order<0){$minus+=1;}
		}
		my $multi=1;
		print "before",join "\t",@asort,"\t";
		if($minus > $plus){@asort=reverse @asort;$multi=-1;}
		#if($ff_order<0){@asort=reverse @asort;}
		print join "\t",@asort,"\n";
		my @all_bsort;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			$asort_order*=$multi;
			my @bsort=sort {$a<=>$b} keys %{$anchor{$asort}};
			if($asort_order<0){@bsort=sort {$b<=>$a} @bsort;}

			for(my $i=0;$i<@bsort;$i+=1){
				my $bsort=$bsort[$i];
				my @tmp=@{$anchor{$asort}{$bsort}};
				push @all_bsort,[@tmp];
				print "tt\t$asort\t$asort_order\t",join "\t",@tmp,"\n";
			}
		}
		my $pre_sp;my @pre_bsort;
		for(my $i=0;$i<@all_bsort;$i+=1){
			my @tmp=@{$all_bsort[$i]};
			my $sp=$tmp[0];
			unless($pre_sp){$pre_sp=$sp;push @pre_bsort,$i;next;}
			if($sp ne $pre_sp){
				#need to output the pre sp record
				my $ast=$pre_bsort[0];
				my @key=sort {$a<=>$b} keys %{$link{$pre_sp}};
				my @atmp=@{$all_bsort[$ast]};
				my $bacsort=$atmp[2];
				my $sp_order=1;
				if($bacsort == $link{$pre_sp}{$key[-1]}[-1]){@key=sort {$b<=>$a} @key;$sp_order=-1;}
				if(@pre_bsort==1 and $ast==0 and $bacsort==$link{$pre_sp}{$key[0]}[1]){
					@key=sort {$b<=>$a} @key;$sp_order=-1;
				}
				foreach my $key(@key){
					my @key1=@{$link{$pre_sp}{$key}};
					my $ctg=shift @key1;
					#if($sp_order>0){@key1=sort {$a<=>$b} @key1;}
					if($sp_order<0){@key1=reverse @key1;}
					#$last_order=$sp_order;
					foreach(@key1){
						print "finalCTG\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$_\ttt\n";
						push @{$bac_sort{$ff}{$bac_sort}},[$ctg,$_];
						$on_bac{$ctg}{$_}=[$pre_sp,];
						$bac_sort+=1;
					}
				}
				@pre_bsort=();
				$pre_sp=$sp;push @pre_bsort,$i;
			}else{
				$pre_sp=$sp;push @pre_bsort,$i;
			}
		}
		#print the last part
		my $ast=$pre_bsort[0];
		my @key=sort {$a<=>$b} keys %{$link{$pre_sp}};
		my @atmp=@{$all_bsort[$ast]};
		my $bacsort=$atmp[2];
		my $sp_order=1;
		if($bacsort == $link{$pre_sp}{$key[-1]}[-1]){@key=sort {$b<=>$a} @key;$sp_order=-1;}
		foreach my $key(@key){
			my @key1=@{$link{$pre_sp}{$key}};
			my $ctg=shift @key1;
			#if($sp_order>0){@key1=sort {$a<=>$b} @key1;}
			if($sp_order<0){@key1=reverse @key1;}
			foreach(@key1){
				print "finalCTG\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$_\ttt\n";
				push @{$bac_sort{$ff}{$bac_sort}},[$ctg,$_];
				$on_bac{$ctg}{$_}=[$pre_sp,];
				$bac_sort+=1;
			}
		}
	
	}
}

#find the missing BAC sort
my %input;my %input_bu;
foreach my $ctg(keys %cpt_rel){
	my %miss;my %out_on;my @on_count;my $max;
	foreach my $ss(sort {$a<=>$b} keys %{$cpt_rel{$ctg}}){
		push @{$out_on{$ss}},$ctg;$max=$ss;
		if($on_bac{$ctg} and $on_bac{$ctg}{$ss}){
			push @{$out_on{$ss}},"on";push @on_count,$ss;
			print "$ctg\t$ss\ton\t$on_bac{$ctg}{$ss}[0]\n";
			next;
		}
		push @{$out_on{$ss}},"out";
		print "$ctg\t$ss\tout\n";
		$miss{$ss}=1;
		#print "$ctg\t$ss\tmiss\n";
	}
	if(@on_count==1 and $lone_ctg_sp{$ctg}){
		#all on one;
		$input_bu{$ctg}{$on_count[0]}=$max;
		print "lone\t$ctg\t$on_count[0]\t$max\n";
	}
	#find a place to locate the miss bac sort and output the last bac sort;some place may need to split by the overlap
	#information because they can not link together
	if($on_bac{$ctg} and keys %miss >0){
		my @key=sort {$a<=>$b} keys %out_on;
		my $st=-1;my $ed=-1;
		my %both;
		for(my $i=0;$i<@key;$i+=1){
			if($out_on{$key[$i]}[1] eq "out"){
				if($st<0){$st=$i;}
			}
			if($out_on{$key[$i]}[1] eq "on" and $st>=0){
				if($ed<0) {
					$ed=$i;
				}
			}
			if($st>=0 and $ed>=0){
				#get out block; start to decide which block it belong to
				my $pre_on=$st-1;if($st==0){$pre_on=-1;}
				my $aft_on=$ed+1;if($ed+1>@key){$aft_on=-1;}
				my $check_gap=0;my $follow;my $follow1;
				if($pre_on<0){
					#no pre
					if($aft_on<0){
						#no anchor
						print "input\t$out_on{$key[$i]}[0]\tun_anchor\n";
					}else{
						#input on aft
						print "input\t$ctg\t$st\t$ed\taft\t$ed\t",join "\t",($st..$ed-1),"\n";
						$input{$ctg}{$ed}=[($st..$ed-1),'pre'];$check_gap=1;$follow=$ed;$follow1=$ed-1;
					}
				}else{
					if($aft_on<0){
						#input on pre
						print "input\t$ctg\t$st\t$ed\tpre\t$st\t",join "\t",($st+1..$ed),"\n";
						$input{$ctg}{$st}=[($st+1..$ed),'aft'];$check_gap=1;$follow=$st;$follow1=$st+1;
					}else{
						#both. check which one is properble. check the on number use the larger one
						print "input\t$ctg\t$st\t$ed\tboth\t";
						#check pre on
						my $p_on=0;
						for(my $j=$st-1;$j>=0;$j-=1){
							if($out_on{$key[$j]}[1] eq 'out'){last;}
							if($out_on{$key[$j]}[1] eq 'on'){$p_on+=1;}
						}
						my $a_on=0;
						for(my $j=$ed;$j<@key;$j+=1){
							if($out_on{$key[$j]}[1] eq 'out'){last;}
							if($out_on{$key[$j]}[1] eq 'on'){$a_on+=1;}
						}
						if($p_on>=$a_on){
							#choose p
							my $pre=$st-1;
							$input{$ctg}{$pre}=[($st..$ed-1),'aft'];$check_gap=1;$follow=$pre;$follow1=$st;
							print "$pre\t",join "\t",($st..$ed-1),"\n";
						}else{
							#choose a
							$input{$ctg}{$ed}=[($st..$ed-1),'pre'];$check_gap=1;$follow=$ed;$follow1=$ed-1;
							print "$ed\t",join "\t",($st..$ed-1),"\n";
						}
					}
				}
				#check gap stat
				if($check_gap){
					my $sp_on=$on_bac{$ctg}{$follow}[0];
					my ($y_or_n,$sp_in_order)=check_gap($follow,$follow1,$sp_on,$ctg,$all_ctg{$ctg},\%all);
					if($y_or_n){push @{$input{$ctg}{$follow}},$sp_in_order;}
					else{
						print "nofollow\t$ctg\t$follow\tcare\t",join "\t",@{$input{$ctg}{$follow}},"\n";
						undef $input{$ctg}{$follow};delete $input{$ctg}{$follow};
					}
				}
				$st=-1;$ed=-1;
			}
			#$st=-1;$ed=-1;
		}
		my $check_gap=0;my $follow;my $follow1;
		if($st>=0 and $ed<0){
			#in last;
			my $pre=$st-1;
			print "input\t$ctg\t$pre\t$key[-1]\tlast\t$pre\t",join "\t",($st..$key[-1]),"\n";
			$input{$ctg}{$pre}=[($st..$key[-1]),'aft'];$check_gap=1;$follow=$pre;$follow1=$st;
		}
		if($check_gap){
			my $sp_on=$on_bac{$ctg}{$follow}[0];
			my ($y_or_n,$sp_in_order)=check_gap($follow,$follow1,$sp_on,$ctg,$all_ctg{$ctg},\%all);
			if($y_or_n){push @{$input{$ctg}{$follow}},$sp_in_order;}
			else{
				print "nofollow\t$ctg\t$follow\tcare\t",join "\t",@{$input{$ctg}{$follow}},"\n";
				undef $input{$ctg}{$follow};delete $input{$ctg}{$follow};
			}
		}
	}
	
}

#output  all  result again
%ex=();
%bac_sort=();%on_bac=();my %last_output;
#foreach my $ff(sort {(keys $final{$b})<=>(keys $final{$a})} keys %final){
foreach my $ff(sort {@{$final{$b}}<=>@{$final{$a}}} keys %final){
	my $bac_sort=1;
	my $last_output=1;
	my @asort=@{$final{$ff}};
	my $ff_order=pop @asort;
	my $out=join "\t",@asort;
	my $po=1;
	foreach my $ii(@asort){
		my ($asort,$order)=split /_/,$ii;
		if($ex{$asort}){$po=0;last;}
		$ex{$asort}=1;
	}
	if($po){
		#this sort will be used
		my $plus=0;my $minus=0;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			if($asort_order>0){$plus+=1;}if($asort_order<0){$minus+=1;}
		}
		my $multi=1;
		if($minus > $plus){@asort=reverse @asort;$multi=-1;}
		#if($ff_order<0){@asort=reverse @asort;}
		my @all_bsort;
		foreach my $ii(@asort){
			my ($asort,$asort_order)=split /_/,$ii;
			$asort_order*=$multi;
			my @bsort=sort {$a<=>$b} keys %{$anchor{$asort}};
			if($asort_order<0){@bsort=sort {$b<=>$a} @bsort;}

			for(my $i=0;$i<@bsort;$i+=1){
				my $bsort=$bsort[$i];
				my @tmp=@{$anchor{$asort}{$bsort}};
				push @all_bsort,[@tmp];
			}
		}
		my $pre_sp;my @pre_bsort;
		for(my $i=0;$i<@all_bsort;$i+=1){
			my @tmp=@{$all_bsort[$i]};
			my $sp=$tmp[0];
			unless($pre_sp){$pre_sp=$sp;push @pre_bsort,$i;next;}
			if($sp ne $pre_sp){
				#need to output the pre sp record
				my $ast=$pre_bsort[0];
				my @key=sort {$a<=>$b} keys %{$link{$pre_sp}};
				my @atmp=@{$all_bsort[$ast]};
				my $bacsort=$atmp[2];
				my $sp_order=1;
				if($bacsort == $link{$pre_sp}{$key[-1]}[-1]){@key=sort {$b<=>$a} @key;$sp_order=-1;}
				if(@pre_bsort==1 and $ast==0 and $bacsort==$link{$pre_sp}{$key[0]}[1]){
					@key=sort {$b<=>$a} @key;$sp_order=-1;
				}
				foreach my $key(@key){
					my @key1=@{$link{$pre_sp}{$key}};
					my $ctg=shift @key1;
					
					#if($sp_order>0){@key1=sort {$a<=>$b} @key1;}
					if($sp_order<0){@key1=reverse @key1;}
					#$last_order=$sp_order;
					if($input_bu{$ctg}){
						#
						my $max=$input_bu{$ctg}{$key1[0]};
						foreach my $tmpss((0..$max)){
							print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$tmpss\tlo\n";
							$last_output{$ff}{$last_output}=[$ctg,$tmpss,$pre_sp];$last_output+=1;
						}
						next;
					}

					foreach my $key1(@key1){
						if($input{$ctg}{$key1}){
							my @pp=@{$input{$ctg}{$key1}};print "test\t$ctg\t",join "\t",@pp,"\n";
							my $pp_sort=1;
							my $spinorder=pop @pp;
							my $last=pop @pp;my $pri;

							if($spinorder*$sp_order<0){@pp=reverse @pp;}
							if($sp_order>0){
								if($last eq 'pre'){
									if($spinorder<0){
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
									}else{
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\t$spinorder\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
									}
								}else{
									if($spinorder<0){
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
									}else{
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
									}
								}
							}else{
								if($last eq 'pre'){
									if($spinorder<0){
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
									}else{
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
									}
								}else{
									if($spinorder<0){
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
									}else{
										foreach my $new(@pp){
											print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
											$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
										}
										print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
										$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
									}
								}
							}
							next;
						}
						print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\tttb\n";
						$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
						push @{$bac_sort{$ff}{$bac_sort}},$ctg,$key1;
						$on_bac{$ctg}{$key1}=[$pre_sp,];
						$bac_sort+=1;
					}
				}
				@pre_bsort=();
				$pre_sp=$sp;push @pre_bsort,$i;
			}else{
				$pre_sp=$sp;push @pre_bsort,$i;
			}
		}
		#print the last part
		my $ast=$pre_bsort[0];
		my @key=sort {$a<=>$b} keys %{$link{$pre_sp}};
		my @atmp=@{$all_bsort[$ast]};
		my $bacsort=$atmp[2];
		my $sp_order=1;
		if($bacsort == $link{$pre_sp}{$key[-1]}[-1]){@key=sort {$b<=>$a} @key;$sp_order=-1;}
		foreach my $key(@key){
			my @key1=@{$link{$pre_sp}{$key}};
			my $ctg=shift @key1;
			#if($sp_order>0){@key1=sort {$a<=>$b} @key1;}
			if($sp_order<0){@key1=reverse @key1;}
			if($input_bu{$ctg}){
				#
				my $max=$input_bu{$ctg}{$key1[0]};
				foreach my $tmpss((0..$max)){
					print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$tmpss\tlo\n";
					$last_output{$ff}{$last_output}=[$ctg,$tmpss,$pre_sp];$last_output+=1;
				}
			}
			foreach my $key1(@key1){
				if($input{$ctg}{$key1}){
					my @pp=@{$input{$ctg}{$key1}};
					my $pp_sort=1;
					my $spinorder=pop @pp;
					my $last=pop @pp;my $pri;
					if($spinorder*$sp_order<0){@pp=reverse @pp;}
					if($sp_order>0){
						if($last eq 'pre'){
							if($spinorder<0){
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
							}else{
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\t$spinorder\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
							}
						}else{
							if($spinorder<0){
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
							}else{
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
							}
						}
					}else{
						if($last eq 'pre'){
							if($spinorder<0){
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
							}else{
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
							}
						}else{
							if($spinorder<0){
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
							}else{
								foreach my $new(@pp){
									print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$new\tin\n";
									$last_output{$ff}{$last_output}=[$ctg,$new,$pre_sp];$last_output+=1;
								}
								print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\ttt\n";
								$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
							}
						}
					}
					next;
				}
				print "haha\t$ff\t$key\t$pre_sp\t$sp_order\t$ctg\t$key1\tttbc\n";
				$last_output{$ff}{$last_output}=[$ctg,$key1,$pre_sp];$last_output+=1;
				push @{$bac_sort{$ff}{$bac_sort}},[$ctg,$key1];
				$on_bac{$ctg}{$key1}=[$pre_sp,];
				$bac_sort+=1;
			}
		}
	}
}

# check CTG split stat
my %block;
my %new_sc;my $cu=1;
foreach my $ff(keys %last_output){
	my @key=sort {$a<=>$b} keys %{$last_output{$ff}};
	my $pre_ctg;
	foreach my $key(@key){
		unless($pre_ctg){$pre_ctg=$last_output{$ff}{$key}[0];}
		if($pre_ctg ne $last_output{$ff}{$key}[0]){$cu+=1;}
		push @{$block{$last_output{$ff}{$key}[0]}{$cu}},$last_output{$ff}{$key}[1];
		$new_sc{$ff}{$cu}{$last_output{$ff}{$key}[0]}=1;
		$pre_ctg=$last_output{$ff}{$key}[0];
	}
}

my %split;
foreach my $ctg(keys %block){
	my @cu=keys %{$block{$ctg}};#block number
	foreach my $cu(@cu){
		my @tmp=@{$block{$ctg}{$cu}};
		print "splitblock\t$cu\t$ctg\t",join "\t",@tmp,"\n";
		my $pre=-1;my $aft;
		for(my $i=1;$i<@tmp;$i+=1){
			my $st=$tmp[$i-1];my $ed=$tmp[$i];
			if(abs($st-$ed)==1){
				if($pre<0){$pre=$st;}
				next;
			}else{
				if($pre<0){$pre=$st;}
				$aft=$st;
				push @{$split{$ctg}{$cu}},[$pre,$aft];
				$pre=-1;
			}
		}
		if($pre>=0){
			push @{$split{$ctg}{$cu}},[$pre,$tmp[-1]];
		}
	}
}

foreach my $ctg(keys %split){
	print "xixi\t$ctg\t";
	foreach my $ccu(sort {$a<=>$b} keys %{$split{$ctg}}){
		print "$ccu\t";
		foreach my $tmp(@{$split{$ctg}{$ccu}}){
			print "$tmp->[0]-$tmp->[1]\t";
		}
	}
	print "\n";
}

#find the split point of the ctg
#%cpt_rel  %ctg_agp
my %take;my %used;
foreach my $ctg(keys %split){
	foreach my $ccu(sort {$a<=>$b} keys %{$split{$ctg}}){
		my $pp=1;my %pp;
		foreach my $tmp(@{$split{$ctg}{$ccu}}){
			my $st=$tmp->[0];my $ed=$tmp->[1];
			$pp{$pp}=[$st,$ed];$pp+=1;
		}
		foreach my $ppp(sort {$pp{$a}[0]<=>$pp{$b}[0]} keys %pp)
		{
			my $st=$pp{$ppp}[0];my $ed=$pp{$ppp}[1];
			($st,$ed)=sort {$a<=>$b} ($st,$ed);
			my %bac;
			foreach(($st..$ed)){
				my $bac=$cpt_rel{$ctg}{$_}[0];
				print "$ctg\t$_\t$bac\n";
				$bac{$bac}=1;
			}
			my $stt=-1;my $edd=0;
			foreach my $ss(sort {$a<=>$b} keys %{$ctg_agp{$ctg}}){
				my @aa=@{$ctg_agp{$ctg}{$ss}};
				if($aa[5]=~/bac/ or $aa[-1]=~/bac/){
					my @bac=grep /bac/,@aa;
					$bac[0]=~s/_\d+$//;
					my $add=0;
					if($aa[-1]=~/bac/){$add=1;}
					print "yaya\t$ctg\t$add\t$bac[0]\t$ss\t";
					if($bac{$bac[0]}){
						print "1\t";
						if($add){
							if($stt<0){$stt=$ss;}
							if($ss>$edd){$edd=$ss;}
						}else{
							if($aa[-1] == 1){
								if($stt<0){$stt=$ss;}
							}
							if($aa[-1] == 2){
								if($ss>$edd){$edd=$ss;}
							}
						}
					}
					print "\n";
				}
			}
			my $flag=0;
			print "$ctg\t$st\t$ed\t$stt\t$edd\n";
			if($stt<0 and $edd==0){
				print STDERR "$ctg\tno assembly\n";
			}
			foreach my $ss(sort {$a<=>$b} keys %{$ctg_agp{$ctg}}){
				if($ss == $stt){
					$flag=1;
					my $out=join "\t",@{$ctg_agp{$ctg}{$ss}};
					$take{$ctg}{$ccu}{$ss}=$out;$used{$ctg}{$ss}=$ccu;
					next;
				}
				if($ss == $edd){
					$flag=0;
					my $out=join "\t",@{$ctg_agp{$ctg}{$ss}};
					$take{$ctg}{$ccu}{$ss}=$out;$used{$ctg}{$ss}=$ccu;
					next;
				}
				if($flag){
					my $out=join "\t",@{$ctg_agp{$ctg}{$ss}};
					$take{$ctg}{$ccu}{$ss}=$out;$used{$ctg}{$ss}=$ccu;
				}
			}
		}
	}
}

#check end. make sure the record after the end is scaffold. and supply the unused record to proper place.
foreach my $ctg(keys %ctg_agp){
	my $pre_ss=-1;my @add;
	foreach my $ss(sort {$a<=>$b} keys %{$ctg_agp{$ctg}}){
		if($used{$ctg}{$ss}){
			$pre_ss=$ss;
			if(@add>0){
				foreach my $ass(@add){
					my $out=join "\t",@{$ctg_agp{$ctg}{$ass}};
					my $ccu=$used{$ctg}{$pre_ss};
					$take{$ctg}{$ccu}{$ass}=$out;
				}
			}
			@add=();
			next;
		}
		my @aa=@{$ctg_agp{$ctg}{$ss}};
		my @bb=grep /bac/,@aa;
		if(@bb==0){
			if(@add>0){push @add,$ss;}
			next;
		}
		else{
			if($ss-$pre_ss>1){
				#add to aft
				push @add,$ss;
			}else{
				#add to this
				my $out=join "\t",@{$ctg_agp{$ctg}{$ss}};
				my $ccu=$used{$ctg}{$pre_ss};
				$take{$ctg}{$ccu}{$ss}=$out;
			}
		}
	}
}


open OUT,">new.agp" or die ;
foreach my $ff(keys %new_sc){
	foreach my $ccu(sort {$a<=>$b} keys %{$new_sc{$ff}}){
		foreach my $ctg(keys %{$new_sc{$ff}{$ccu}}){
			foreach my $re(sort {$a<=>$b} keys %{$take{$ctg}{$ccu}}){
				print OUT "new$ff\t$ccu\t$take{$ctg}{$ccu}{$re}\n";
			}
			print OUT "new$ff\t$ccu\t$ctg\t1\t100\t1\tN\t100\tscaffold\tno\tna\n";
		}
		print OUT "new$ff\t$ccu\tadd\t1\t100\t1\tN\t100\tscaffold\tno\tna\n";

	}
}
close OUT;





#sub routine
sub check_gap
{
	my ($bac_sort,$bac_sort_fow,$sp,$ctg,$agp_hash,$all_hash)=@_;
	my @key=sort {$a<=>$b} keys %{$all_hash->{$sp}};
	my $pos;
	for(my $i=0;$i<@key;$i++){
		if(@{$all_hash->{$sp}{$key[$i]}}>10 and $all_hash->{$sp}{$key[$i]}[9] eq $ctg and $all_hash->{$sp}{$key[$i]}[10] == $bac_sort){
			$pos=$i;last;
		}
	}
	my $pre=0;my $aft=0;
	if($pos-2>=0 and $all_hash->{$sp}{$key[$pos-2]} and $all_hash->{$sp}{$key[$pos-2]}[9] eq $ctg){
		$pre=1;
	}
	if($pos+2<@key and $all_hash->{$sp}{$key[$pos+2]} and $all_hash->{$sp}{$key[$pos+2]}[9] eq $ctg){
		$aft=1;
	}
	if($pre){
		#not pre;
		my $gap_length;
		my $spinorder=$bac_sort-$all_hash->{$sp}{$key[$pos-2]}[10];
		$spinorder=$spinorder/abs($spinorder);
		#print "$spinorder\n";
		if($pos+1<@key and $all_hash->{$sp}{$key[$pos+1]}){
			$gap_length=$all_hash->{$sp}{$key[$pos+1]}[5];
			if($gap_length<10000){
				return -1,$spinorder;
			}else{
				return 1,$spinorder;
			}
		}else{
			#in end
			return 1,$spinorder;
		}
		
	}elsif($aft){
		#not aft
		my $gap_length;
		print STDERR "$sp\t$ctg\t$bac_sort\t$all_hash->{$sp}{$key[$pos+2]}[10]\t$pos\t@key\t$aft\t";
		my $spinorder=$all_hash->{$sp}{$key[$pos+2]}[10]-$bac_sort;$spinorder=$spinorder/abs($spinorder);
		print STDERR "$spinorder\n";
		if($pos-1>=0 and $all_hash->{$sp}{$key[$pos-1]}){
			$gap_length=$all_hash->{$sp}{$key[$pos-1]}[5];
			if($gap_length<10000){
				return -1,$spinorder;
			}else{
				return 1,$spinorder;
			}
		}else{
			#in end
			return 1,$spinorder;
		}
	}else{
		#single
		my $gap_length_aft=0;
		if($pos+1<@key and $all_hash->{$sp}{$key[$pos+1]}){
			$gap_length_aft=$all_hash->{$sp}{$key[$pos+1]}[5];
		}
		my $gap_length_pre=0;
		if($pos-1>=0 and $all_hash->{$sp}{$key[$pos-1]}){
			$gap_length_pre=$all_hash->{$sp}{$key[$pos-1]}[5];
		}
		if($gap_length_pre and $gap_length_aft){
			if($gap_length_aft >=$gap_length_pre){
				return 1,1;
			}else{
				return 1,-1;
			}
		}else{
			if($gap_length_pre and $gap_length_pre>10000){
				if($bac_sort>$bac_sort_fow){
					return 1,1;
				}else{
					return 1,-1;
				}
			}
			if($gap_length_aft and $gap_length_aft>10000){
				if($bac_sort>$bac_sort_fow){
					return 1,-1;
				}else{
					return 1,1;
				}
			}
			return -1,1;
		}
	}
	

}




sub Find_Path
{
	my ($type,$order,$hash,$mer_hash,$asort,$bsort,$use)=@_;
	my $pre_sp=$hash->{$asort}{$bsort}[0];
	#print "start $pre_sp\t$asort\t$bsort\n";
	my $nextasort;my $nextbsort;
	my @back;
	$use->{$asort}=1;
	my %dis;
	my $raw_order=$order;my $next_order;
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
			if($type==1){$next_order=-1;}
			elsif($type==-1){$next_order=1;}
		}else{
			$nextb=$next[0];
			if($type==1){$next_order=1;}
			elsif($type==-1){$next_order=-1;}
		}
		if(@next==1){
			$next_order=1;$use->{$nextasort}+=1;
			#same sp;return info to discard this anchor
			$dis{$nextasort}=1;
		}
		if($use->{$nextasort}){
			;
		}else{
			#$tht=$hash->{$nextasort}{$nextb}[3];
			@back=Find_Path($type,$next_order,$hash,$mer_hash,$nextasort,$nextb,$use);
		}
	}
	my $out="$asort";if($raw_order){$out.="_$raw_order";}
	push @back,$out;
	if(keys %dis>0){
		#push @back,"dis";
		foreach(keys %dis){
			#push @back,$_;
		}
	}
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


