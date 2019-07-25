//###################################################################
//### Name:        ShowerProduedPtrsHolder                        ### 
//### Author:      Dominic Barker                                 ###
//### Date:        15.07.19                                       ###
//### Description: Class to holder the unique ptrs required to    ###
//###              produce data products. Used in SBNShower_modle ###
//###              and corresponding tools.                       ###
//###################################################################

#ifndef ShowerProduedPtrsHolder_HH
#define ShowerProduedPtrsHolder_HH

//Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "larreco/ShowerFinder/ShowerElementHolder.hh"
#include "larreco/ShowerFinder/ShowerProduedPtrsHolder.hh"

//C++ Includes 
#include <iostream> 
#include <map>
#include <string> 
#include <memory>

namespace reco {
  namespace shower {
    class ShowerUniqueProduerPtrBase;
    template <class T> class ShowerUniqueProductPtr;
    template <class T> class ShowerUniqueAssnPtr;
    class ShowerPtrMakerBase;
    template <class T> class ShowerPtrMaker;
    class ShowerProduedPtrsHolder;
  }
}

//Templates to check if a type is a std::vector
template <class N>
struct is_vector { static const int value = 0; };

template <class N, class A>
struct is_vector<std::vector<N, A> > { static const int value = 1; };

//Template to check if its an assn 
template <class N>
struct is_assn { static const int value = 0; }; 

template  <class L, class R, class D>    
struct is_assn<art::Assns< L, R, D > >{ static const int value = 1; };    

//Template to carry the type becuase functions can't, annoyingly, be partially specialised
template <typename T>
struct type{};


//Base class for the class that holds the ptr.
class reco::shower::ShowerUniqueProduerPtrBase {

public: 

  virtual ~ShowerUniqueProduerPtrBase() noexcept = default;

  virtual void reset() = 0;

  virtual void AddDataProduct(reco::shower::ShowerElementHolder& selement_holder, std::string Name) = 0;

  virtual void MoveToEvent(art::Event& evt, std::string Name) = 0;

};

//Class that holds a unique ptr for the product. This is what is stored in the map. The product is put into
//the event as a vector so the this holder maintains this unique ptr and the actions to manipulate it.
template <class T>
class reco::shower::ShowerUniqueProductPtr<std::vector<T> >: public reco::shower::ShowerUniqueProduerPtrBase{

public:

  ShowerUniqueProductPtr<std::vector<T> >(){
    ptr = 1;
    showeruniqueptr = std::make_unique<std::vector<T> >();
  }

  //Get the unique ptr for the data product.
  std::unique_ptr<T>& GetPtr(){
    if(ptr){
      return showeruniqueptr;
    }
    else{
      throw cet::exception("ShowerUniqueProduerPtr") << "Element does not exist" << std::endl;
      return showeruniqueptr;
    }
  }

  void reset() override {
    showeruniqueptr.reset(new std::vector<T>());
  }

  //Add a data product on to the vector that will be added to the event.
  void AddDataProduct(reco::shower::ShowerElementHolder& selement_holder, std::string Name) override {
    T product;
    selement_holder.GetElement(Name, product);
    showeruniqueptr->push_back(product);
  }

  //Final thing to do move to the event.
  void MoveToEvent(art::Event& evt, std::string Name) override {
    evt.put(std::move(showeruniqueptr),Name);
  }

private:
  std::unique_ptr<std::vector<T> > showeruniqueptr;
  int ptr;
};


//Class that holds a unique ptr for the association. This is what is stored in the map. The association is put into
//the event as a vector so the this holder maintains this unique ptr and the actions to manipulate it.
//I guess if the product if a product is unique to the event then this holder will deal with it.
//I have tried to be smart and I don't think it not being an association is a problem as 
//long as the user is smart. 
template <class T>
class reco::shower::ShowerUniqueAssnPtr: public reco::shower::ShowerUniqueProduerPtrBase{

public:

  ShowerUniqueAssnPtr(){
    ptr = 1;
    showeruniqueptr = std::make_unique<T>();
  }

  //Get the ptr to the association.
  std::unique_ptr<T>& GetPtr(){
    if(ptr){
      return showeruniqueptr;
    }
    else{
      throw cet::exception("ShowerUniqueAssnPtr") << "Element does not exist" << std::endl;
      return showeruniqueptr;
    }
  }

  void reset() override {
    showeruniqueptr.reset(new T());
  }

  //place the association to the event.
  void MoveToEvent(art::Event& evt, std::string Name) override {
    evt.put(std::move(showeruniqueptr), Name);
  }

  //Not need but the compiler complains if its not here.
  void AddDataProduct(reco::shower::ShowerElementHolder& selement_holder, std::string Name) override {std::cout << "should not be here" << std::endl;}


private:
  std::unique_ptr<T> showeruniqueptr;
  int ptr;
};


//Base class to hold the pointer makers. This interacts the with the module a little differently 
//as the ptr makers do not work within the tools. The holds and the ptrmaker and provides set and 
//get functions so that the usr can access art::Ptrs for the products of the module and associate 
//them to other things.
class reco::shower::ShowerPtrMakerBase{

public:

  virtual ~ShowerPtrMakerBase() noexcept = default;
  
  virtual bool CheckPtrMaker() = 0;
    
  virtual void SetPtrMaker(art::Event& evt, std::string InstanceName) = 0;

  virtual void Reset() = 0;

};

//Derived class - See avove
template<class T> 
class reco::shower::ShowerPtrMaker : public ShowerPtrMakerBase{

public:

  ShowerPtrMaker(art::Event& evt, std::string InstanceName){
    ptrmaker = new std::unique_ptr<art::PtrMaker<T> >(new art::PtrMaker<T>(evt,InstanceName));
    ptr = 1;
  }

  ShowerPtrMaker(){
    ptrmaker = nullptr;
    ptr = 0;
  }

  bool CheckPtrMaker() override {
    if(ptr){
      return true;
    }
    return false;
  }

  //Return the ptr maker. Probably never needed.
  art::PtrMaker<T>& GetPtrMaker(){
    if(ptr){
      if(ptrmaker == NULL){
	throw cet::exception("ShowerPtrMaker") << "Ptr maker ptr is null" << std::endl;
      }
      return *ptrmaker;
    }
    throw cet::exception("ShowerPtrMaker") << "Trying to get a  ptrmaker that does not exists" << std::endl;
    return *ptrmaker;
  }

  //Return the art ptr that the module produces corresponding the index iter
  art::Ptr<T> GetArtPtr(int iter){
    if(ptr){
      if(ptrmaker == NULL){
	throw cet::exception("ShowerPtrMaker") << "Ptr maker ptr is null" << std::endl;
      }
      return (*ptrmaker)(iter);
    }
    throw cet::exception("ShowerPtrMaker") << "Trying to get a  ptrmaker that does not exists" << std::endl;
    return (*ptrmaker)(0);
  }

  //Set the ptr maker this is reset at the start of the event. 
  void SetPtrMaker(art::Event& evt, std::string InstanceName) override {
    ptrmaker.reset(new art::PtrMaker<T>(evt,InstanceName));
    ptr = 1;
  }

  void Reset() override {
    if(!ptr){
      throw cet::exception("ShowerPtrMaker") << "Trying to reset ptr but it has not been set in the first place" << std::endl;
    }
    ptrmaker.reset(nullptr);
    ptr = 0;
  }

private:

  std::unique_ptr<art::PtrMaker<T> > ptrmaker; 
  int ptr;
};

//Class that holds all the unique ptrs and the ptr makers. It is what the tools and module use 
//to access the above class elements. The end user case will see the user not interact with this 
// directly. 
class reco::shower::ShowerProduedPtrsHolder {

public: 
  
  //Initialise the a unique ptr in the map. This will be added to the event.
  template <class T>
  int SetShowerUniqueProduerPtr(std::string InstanceName, type<T>){
      //Add to the assns
      if(showerassnPtrs.find(InstanceName) != showerassnPtrs.end()){
	mf::LogError("ShowerProduedPtrsHolder") << "Trying to set ptr: " << InstanceName << ". Instance already exists" << std::endl;
	return 1;
      }
      showerassnPtrs[InstanceName] = std::unique_ptr<reco::shower::ShowerUniqueAssnPtr<T> >(new reco::shower::ShowerUniqueAssnPtr<T>());
      return 0;
  }
  
  //Set the unique ptr. The unique ptr will be filled into the event.
  template <class T>
  int SetShowerUniqueProduerPtr(std::string InstanceName, type<std::vector<T> >){

    //Then add the products
    if(showerproductPtrs.find(InstanceName) != showerproductPtrs.end()){
      mf::LogError("ShowerProduedPtrsHolder") << "Trying to set ptr: " << InstanceName << ". Instance already exists" << std::endl;
      return 1;
    }
    if(showerPtrMakers.find(InstanceName) != showerPtrMakers.end()){
      throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker already exist. It should not be set again" << std::endl;
    }
    showerPtrMakers[InstanceName]   = std::unique_ptr<reco::shower::ShowerPtrMaker<T> >(new reco::shower::ShowerPtrMaker<T>());
    showerproductPtrs[InstanceName] = std::unique_ptr<reco::shower::ShowerUniqueProductPtr<std::vector<T > > >(new reco::shower::ShowerUniqueProductPtr<std::vector<T> >());
    return 0;
  }

  
  //Checks if the ptr exists
  bool CheckUniqueProduerPtr(std::string InstanceName){
    if(showerproductPtrs.find(InstanceName) != showerproductPtrs.end()){
      return true;
    }
    if(showerassnPtrs.find(InstanceName) != showerassnPtrs.end()){
      return true;
    }
    return false;
  }

  //Reset the ptrs;
  void reset(){
    for(auto const& showerptr: showerproductPtrs){
      (showerptr.second)->reset();
    }
    for(auto const& showerptr: showerassnPtrs){
      (showerptr.second)->reset();
    }
  }
 
  //Add any data products that are produced by the module to the unique ptr it corresponds to
  //This is done by matching strings in the element holder and the ptr holder. Hence these
  //must match. This is a global command done in the module.
  void AddDataProducts(reco::shower::ShowerElementHolder& selement_holder){
    for(auto const& showerproductPtr: showerproductPtrs){
      (showerproductPtr.second)->AddDataProduct(selement_holder, showerproductPtr.first);
    }
  }

  //Global command to move all products into the event. This is done in the module.
  void MoveAllToEvent(art::Event& evt){
    for(auto const& showerproductPtr: showerproductPtrs){
      (showerproductPtr.second)->MoveToEvent(evt,showerproductPtr.first);
    }
    for(auto const& showerassnPtr: showerassnPtrs){
      (showerassnPtr.second)->MoveToEvent(evt,showerassnPtr.first);
    }
  }
  
  //This returns the unique ptr. This is a legacy code.
  template <class T>  
  T& GetPtr(std::string InstanceName){
    if(showerproductPtrs.find(InstanceName) != showerproductPtrs.end()){ 
      reco::shower::ShowerUniqueProductPtr<T>* prod = dynamic_cast<reco::shower::ShowerUniqueProductPtr<T> *>(showerproductPtrs[InstanceName].get());
      return prod->GetPtr();
    }
    else if(showerassnPtrs.find(InstanceName) != showerassnPtrs.end()){ 
      reco::shower::ShowerUniqueAssnPtr<T>* assn = dynamic_cast<reco::shower::ShowerUniqueAssnPtr<T> *>(showerassnPtrs[InstanceName].get());
      return assn->GetPtr();   
    }
    else{
      throw cet::exception("ShowerProduedPtrsHolder") << "Element does not exist" << std::endl;
      return T(); 
    }
  }
  
  //Wrapper so that the use the addSingle command for the association. Add A and B to the association just 
  //as if add single add.
  template <class T, class A, class B>
  void AddSingle(A& a, B& b, std::string InstanceName){
    if(showerassnPtrs.find(InstanceName) == showerassnPtrs.end()){ 
      throw cet::exception("ShowerProduedPtrsHolder") << "Element does not exist" << std::endl;
      return;
    }
    if(!is_assn<T>::value){
      throw cet::exception("ShowerProduedPtrsHolder") << "Element is not an assoication" << std::endl;
      return;
    }
    reco::shower::ShowerUniqueAssnPtr<T>* assnptr = dynamic_cast<reco::shower::ShowerUniqueAssnPtr<T> *>(showerassnPtrs[InstanceName].get());
    
    T* assn = dynamic_cast<T*>(assnptr->GetPtr().get()); 
    assn->addSingle(a,b);
    return;
  }  

  //Initialise the ptr makers. This is done at the the start of the module.
  void SetPtrMakers(art::Event& evt){
    for(auto const&  showerPtrMaker:  showerPtrMakers){ 
      if(showerPtrMakers.find(showerPtrMaker.first) == showerPtrMakers.end()){
	throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker was empty. This is concerning" << std::endl;
      }
      showerPtrMakers[showerPtrMaker.first]->SetPtrMaker(evt,showerPtrMaker.first);
    }
  }

  //Wrapper to access a particle PtrMaker. This is legacy as is not used.
  template <class T> 
  art::PtrMaker<T>& GetPtrMaker(std::string InstanceName){
    if(showerPtrMakers.find(InstanceName) == showerPtrMakers.end()){
      throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker does not exist" << std::endl;
    }
    else{
      if(!showerPtrMakers[InstanceName]->CheckPtrMaker()){
	throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker is not set" << std::endl;
      }
      reco::shower::ShowerPtrMaker<T>* ptrmaker = dynamic_cast<reco::shower::ShowerPtrMaker<T> *>(showerassnPtrs[InstanceName].get());
      return ptrmaker->GetPtrMaker();
    }
  }

  //Wrapper to return to the the user the art ptr corresponding the index iter. 
  template <class T> 
  art::Ptr<T> GetArtPtr(std::string InstanceName, int iter){
    if(showerPtrMakers.find(InstanceName) == showerPtrMakers.end()){
      throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker does not exist" << std::endl;
    }
    else{
      if(!showerPtrMakers[InstanceName]->CheckPtrMaker()){
	throw cet::exception("ShowerProduedPtrsHolder") << "PtrMaker is not set" << std::endl;
      }
      reco::shower::ShowerPtrMaker<T>* ptrmaker = dynamic_cast<reco::shower::ShowerPtrMaker<T> *>(showerPtrMakers[InstanceName].get());
      return ptrmaker->GetArtPtr(iter);
    }
  }

  //Legacy not used.
  void ResetPtrMakers(){
    for(auto const& showerPtrMaker: showerPtrMakers){
      (showerPtrMaker.second)->Reset();
    }
  }


private:
  std::map<std::string,std::unique_ptr<reco::shower::ShowerUniqueProduerPtrBase > > showerproductPtrs;
  std::map<std::string,std::unique_ptr<reco::shower::ShowerUniqueProduerPtrBase > > showerassnPtrs;
  std::map<std::string,std::unique_ptr<reco::shower::ShowerPtrMakerBase> > showerPtrMakers;
};



#endif
