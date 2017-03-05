#ifndef _Chain_Defined_
#define _Chain_Defined_
namespace FractalSpace
{
  class Chain
  {
    bool chain_started;
    bool chain_finished;
    bool ok_start_chain;
  public:
    list <Group*> list_groups;
    static int number_chains;
    Chain()
    {
      assert(this);
      chain_started=false;
      chain_finished=false;
      ok_start_chain=false;
      number_chains++;
      //    cout << "Making Chain " << this  << " " << number_chains << endl;
    }
    ~Chain()
    {    
      number_chains--;
      //    cout << "Ending Chain " << this << " " << number_chains << endl;
    }
    void set_ok_start_chain(const bool& haha)
    {
      ok_start_chain=haha;
    }
    bool get_ok_start_chain()
    {
      return ok_start_chain;
    }
    bool get_chain_finished()
    {
      return chain_finished;
    }
    void set_chain_finished(const bool& f)
    {
      chain_finished=f;
    }
    bool get_chain_started()
    {
      return chain_started;
    }
    void set_chain_started(const bool& f)
    {
      chain_started=f;
    }
    static void show_chains(list <Chain*>& list_chains)
    {
      for(list <Chain*>::const_iterator chain_itr=list_chains.begin();chain_itr != list_chains.end();++chain_itr)
	{
	  Chain& chain=**chain_itr;
	  cout << " showing chain " << &chain << endl;
	  for(list <Group*>::const_iterator group_itr=chain.list_groups.begin();group_itr != chain.list_groups.end();++group_itr)
	    {
	      Group& group=**group_itr;
	      cout << "showing group " << &group << " " << group.get_level() << " " << group.get_p_mother_group() << endl;
	    }
	}
    }
  }
    ;
}
#endif
