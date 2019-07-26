//###################################################################
//### Name:        ShowerElementHolder                           ### 
//### Author:      Dominic Barker                                 ###
//### Date:        15.07.19                                       ###
//### Description: Class to holder the standard shower property   ###
//###              information. Used in SBNShower_modle and       ###
//###              corresponding tools.                           ###
//###################################################################

#ifndef ShowerElementHolder_HH
#define ShowerElementHolder_HH

//Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

//C++ Inlcudes 
#include <iostream> 
#include <map>
#include <string> 
#include <memory>

namespace reco {
  namespace shower {
    class ShowerElementBase;
    template <class T> class ShowerElementAccessor;     
    template <class T> class ShowerDataProduct; 
    template <class T, class T2> class ShowerProperty;
    class ShowerElementHolder;
  }
}

class reco::shower::ShowerElementBase {


public:

  virtual ~ShowerElementBase() noexcept = default;
  
  virtual bool CheckIfSet(){
    throw cet::exception("ShowerElementHolder") << "Trying to check an element that is not a product" << std::endl;
  }
  virtual void SetCheckIfSet(bool check){
    throw cet::exception("ShowerElementHolder") << "Trying to set an element that is not a product" << std::endl;
  }

  bool CheckShowerElement(){
    if(elementPtr) return true;
    else return false;
  }

  void Clear(){
    elementPtr    = 0;
  }


protected:

  int elementPtr;
  
};

//This is a template class which holds a shower property. This holds any object e.g. std::vector<double>, double, TVector3 
//and holds various information which helps the showerproperty holder access the elements. A user should not require any part 
//of this class.  
template <class T>
class reco::shower::ShowerElementAccessor : public reco::shower::ShowerElementBase {

public: 

  void SetShowerElement(T& Element){
    element = Element;
    this->elementPtr = 1;
  }
  
  int GetShowerElement(T& Element){
    if(this->elementPtr){
      Element = element;
      return 0;
    }
    else{
      return 1;
    }
  }

  T GetShowerElement(){
    if(this->elementPtr){
      return element;
    }
    else{
      return T();
    }
  } 
    
protected:
  T   element; 
};

//This class holds shower data products which have the potential to be saved in the art::Event e.g. recob::Track. Note the product itself must be store in the element holder as the object will be destoryed in the CalculateProperty Section otherwise. Associtations can be made during Calculate property tool stage. 
template <class T> 
class reco::shower::ShowerDataProduct : public reco::shower::ShowerElementAccessor<T>{

public:
  
  ShowerDataProduct(){
    this->elementPtr      = 1;
    checkifset            = false;
  }

  ShowerDataProduct(T& Element, bool Checkifset){
    this->elementPtr      = 1;
    this->element         = Element;
    checkifset            =  Checkifset;
  }


  void Clear(){
    this->element = T();
    this->elementPtr    = 0;
  }

  bool CheckIfSet(){
    return checkifset;
  }

  void SetCheckIfSet(bool Checkifset){
    checkifset = Checkifset;
  }

  private:
  bool checkifset;
};

//This class holds shower properties e.g. ShowerDirection. The user must define the associated error  
template <class T, class T2>
class reco::shower::ShowerProperty : public reco::shower::ShowerElementAccessor<T>{

public:

  ShowerProperty(){
    this->elementPtr = 1;
  }
  
  ShowerProperty(T& Element, T2& ElementErr){
    this->elementPtr = 1;
    this->element    = Element;
    propertyErr      = ElementErr;
  }
  

  int GetShowerPropertyError(T2& ElementErr){
    if(this->elementPtr){
      ElementErr = propertyErr;
      return 0;
    }
    else{
      return 1;
    }
  }

  void SetShowerProperty(T& Element, T2& ElementErr){
    this->element    = Element;
    this->elementPtr = 1;
    propertyErr      = ElementErr;
  }

  void Clear(){
    this->element = T();
    this->elementPtr = 0;
  }
  
private:
  T2   propertyErr;

};


//Class to holder all the reco::shower::ShowerElement objects. This is essentially a map from a string the object so people can 
//add an object in a tool and get it back later. 
class reco::shower::ShowerElementHolder{

public:
  
  //Getter function for accessing the shower property e..g the direction ShowerElementHolder.GetElement("MyShowerValue"); The name is used access the value and precise names are required for a complete shower in sbnshower: ShowerStartPosition, ShowerDirection, ShowerEnergy ,ShowerdEdx.
  template <class T >
  int GetElement(std::string Name, T& Element){
    if(showerproperties.find(Name) != showerproperties.end()){
      if(showerproperties[Name]->CheckShowerElement()){
	  reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerproperties[Name].get());
	  showerprop->GetShowerElement(Element);
	  return 0;
	}
      else{
	mf::LogWarning("ShowerElementHolder") << "Trying to get Element" << Name << ". This elment has not been filled" << std::endl;
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerproperties[Name].get());
	showerprop->GetShowerElement(Element);
	return 1;
      }
    }
    else if(showerdataproducts.find(Name) != showerdataproducts.end()){
      if(showerdataproducts[Name]->CheckShowerElement()){
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerdataproducts[Name].get());
	showerprop->GetShowerElement(Element);
	return 0;
      }
      else{
	mf::LogWarning("ShowerElementHolder") << "Trying to get Element" << Name << ". This elment has not been filled" << std::endl;
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerdataproducts[Name].get());
	showerprop->GetShowerElement(Element);
	return 1;
      }
    }
    throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return 1;
  }

  //Alternative get function that returns the object. Not recommended.
  template <class T > 
  T GetElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){ 
      if(showerproperties[Name]->CheckShowerElement()){
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerproperties[Name].get());
	return showerprop->GetShowerElement();
      }
    }
    else if(showerdataproducts.find(Name) != showerdataproducts.end()){   
      if(showerdataproducts[Name]->CheckShowerElement()){ 
	reco::shower::ShowerElementAccessor<T> *showerprop = dynamic_cast<reco::shower::ShowerElementAccessor<T> *>(showerdataproducts[Name].get());
	return showerprop->GetShowerElement();
      }
    }	
    throw cet::exception("ShowerElementHolder") << "Trying to get Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return 1;
  }
  
  //Getter function for accessing the shower property error e.g the direction ShowerElementHolder.GetElement("MyShowerValue");
  template <class T, class T2>
  int GetElementAndError(std::string Name, T& Element,  T2& ElementErr){
    if(showerproperties.find(Name) == showerproperties.end()){
      mf::LogError("ShowerElementHolder") << "Trying to get Element Errorr: " << Name << ". This elment does not exist in the element holder" << std::endl;
      return 1;
    }
    reco::shower::ShowerProperty<T,T2> *showerprop = dynamic_cast<reco::shower::ShowerProperty<T,T2> *>(showerproperties[Name].get());
    showerprop->GetShowerElement(Element);  
    showerprop->GetShowerPropertyError(ElementErr);
    return 0;
  }


  //This sets the value of the data product. Just give a name and a object
  //e.g. TVector3 ShowerElementHolder.SetElement((TVector3) StartPosition, "StartPosition");
  template <class T>
  void SetElement(T& dataproduct, std::string Name, bool checkifset=false){
    
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      reco::shower::ShowerDataProduct<T>* showerdataprod = dynamic_cast<reco::shower::ShowerDataProduct<T> *>(showerdataproducts[Name].get());
      showerdataprod->SetShowerElement(dataproduct);
      showerdataprod->SetCheckIfSet(checkifset);
      return;
    }
    else{
      showerdataproducts[Name] = std::unique_ptr<reco::shower::ShowerDataProduct<T> >(new reco::shower::ShowerDataProduct<T>(dataproduct,checkifset));
      return;
    }
  }

  //This sets the value of the property. Just give a name and a object
  //e.g. TVector3 ShowerElementHolder.SetElement((art::Ptr<recob::Track>) track, "StartPosition", save);
  template <class T, class T2>
  void SetElement(T& propertyval, T2& propertyvalerror, std::string Name){

    if(showerproperties.find(Name) != showerproperties.end()){
      reco::shower::ShowerProperty<T,T2>* showerprop = dynamic_cast<reco::shower::ShowerProperty<T,T2> *>(showerproperties[Name].get());
      showerprop->SetShowerProperty(propertyval,propertyvalerror);
    return;
    }
    else{
      showerproperties[Name] = std::unique_ptr<reco::shower::ShowerProperty<T,T2> >(new reco::shower::ShowerProperty<T,T2>(propertyval,propertyvalerror));
      return;
    }
  }

  //Check that a property is filled 
  bool CheckElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){
      return showerproperties[Name]->CheckShowerElement();
    }
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      return showerdataproducts[Name]->CheckShowerElement();
    }
    return false;
  }

  //Check All the properties
  bool CheckAllElements(){
    bool checked = true;
    for(auto const& showerprop: showerproperties){
     checked *= showerprop.second->CheckShowerElement();
    }
    for(auto const& showerdataprod: showerdataproducts){
      checked *= showerdataprod.second->CheckShowerElement();
    }
    return checked;
  }

  
  //Clear Fucntion. This does not delete the element.
  void ClearElement(std::string Name){
    if(showerproperties.find(Name) != showerproperties.end()){
      return showerproperties[Name]->Clear();
    }
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      return showerdataproducts[Name]->Clear();
    }
    mf::LogError("ShowerElementHolder") << "Trying to clear Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //Clear all the shower properties. This does not delete the element.
  void ClearAll(){
    for(auto const& showerprop: showerproperties){
      (showerprop.second)->Clear();
    }
    for(auto const& showerdataproducts: showerdataproducts){
      (showerdataproducts.second)->Clear();
    }
  }

  //Find if the product is one what is being stored.
  bool CheckElementSavTag(std::string Name){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
    return showerdataproducts[Name]->CheckIfSet();
    }
    return false;
  }
 
  //Delete a product. I see no reason for it.
  void DeleteElement(std::string Name){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      showerdataproducts[Name].reset(nullptr);
      return;
    }
    if(showerproperties.find(Name) != showerproperties.end()){
      showerproperties[Name].reset(nullptr);
      return;
    }
    mf::LogError("ShowerElementHolder") << "Trying to delete Element: " << Name << ". This element does not exist in the element holder" << std::endl;
    return;
  }

  //Set the indicator saying if the shower is going to be stored.
  void SetElementSaveTag(std::string Name, bool checkelement){
    if(showerdataproducts.find(Name) != showerdataproducts.end()){
      showerdataproducts[Name]->SetCheckIfSet(checkelement);
      return;
    }
    mf::LogError("ShowerElementHolder") << "Trying set the checking of the data product: " << Name << ". This data product does not exist in the element holder" << std::endl;
    return;
  }

  //Set the shower number. This is required the association making.
  void SetShowerNumber(int shower_iter){
    showernumber = shower_iter;
  }
  
  //Get the shower number.
  int GetShowerNumber(){
    return showernumber;
  }

private:

  //Storage for all the shower properties.
  std::map<std::string,std::unique_ptr<reco::shower::ShowerElementBase> > showerproperties;
  
  //Storage for all the data products
  std::map<std::string,std::unique_ptr<reco::shower::ShowerElementBase> > showerdataproducts;

  //Shower ID number. Use this to set ptr makers.
  int showernumber;

};


#endif
